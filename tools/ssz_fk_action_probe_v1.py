#!/usr/bin/env python3
"""SSZ-FK action-consistency research probe (separate from frozen Run 11).

A radial point-particle + locally measured elastic and pinning energies on a
fixed static metric.  No claim about a microscopic neutron-star substrate.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment, minimize


@dataclass(frozen=True)
class Metric:
    r_s: float = 4.4
    origin: float = 12.0
    flat: bool = False

    @property
    def A(self) -> float:
        return 0.0 if self.flat else 0.5 * self.r_s

    def D(self, r):
        r = np.asarray(r)
        return r / (r + self.A)

    def Dp(self, r):
        r = np.asarray(r)
        return self.A / (r + self.A) ** 2

    def Dpp(self, r):
        r = np.asarray(r)
        return -2 * self.A / (r + self.A) ** 3

    def s(self, r):
        r = np.asarray(r)
        return 1 + self.A / r

    def sp(self, r):
        r = np.asarray(r)
        return -self.A / r**2

    def ell(self, r):
        r = np.asarray(r)
        return r - self.origin + self.A * np.log(r / self.origin)


@dataclass(frozen=True)
class Parameters:
    sites: int = 40
    R: float = 12.0
    dr: float = 0.05
    r_s: float = 4.4
    kappa: float = 25.0
    m: float = 1.0
    V0: float = 0.0003
    c: float = 1.0
    delta: float = (np.sqrt(5.0) - 2.0) / 4.0
    natural_bonds: str = "background"  # or "constant"
    include_tidal: bool = True
    phase_mode: str = "proper_integral"  # literal_sr | frozen_local
    mass_mode: str = "relativistic_action"  # legacy_m_over_d2


class SSZFKAction:
    """Radial coordinate u_i; r_i = R + i*dr + u_i in *one* coordinate.

    Low-v action, coordinate time:
      L = 1/2 sum m*s(r_i)^2/D(r_i) * dot(r_i)^2 - U
      U = sum_i m*c^2*[D(r_i)-D(q_i)-D'(q_i)*(r_i-q_i)]
        + sum_b D(r_b_mid)*kappa/2*[ell(r_{i+1})-ell(r_i)-ell0_i]^2
        + sum_i D(r_i)*V0*[1-cos(Q*ell(r_i))]
    The linear support potential cancels the background static gravity force.
    These are explicit modeling choices, not a derived relativistic crust EOS.
    """
    def __init__(self, params: Parameters, flat=False):
        self.p = params
        if params.sites < 3 or params.dr <= 0 or params.m <= 0:
            raise ValueError("Invalid parameters")
        if params.phase_mode not in {"proper_integral", "literal_sr", "frozen_local"}:
            raise ValueError("Unknown phase prescription")
        if params.mass_mode not in {"relativistic_action", "legacy_m_over_d2"}:
            raise ValueError("Unknown kinetic prescription")
        if params.natural_bonds not in {"background", "constant"}:
            raise ValueError("Unknown natural-bond prescription")
        self.g = Metric(params.r_s, params.R, flat=flat)
        self.q = params.R + params.dr * np.arange(params.sites)
        self.Q = 2 * np.pi / (params.dr * (1 + params.delta))
        self.l0 = (np.diff(self.g.ell(self.q)) if params.natural_bonds == "background"
                   else np.full(params.sites - 1, params.dr))

    def position(self, u):
        x = self.q + np.asarray(u)
        if np.min(x) <= 0 or np.min(np.diff(x)) <= 0:
            raise ValueError("Sites must remain radially ordered and positive")
        return x

    def mass(self, u):
        x = self.position(u)
        if self.p.mass_mode == "legacy_m_over_d2":
            return self.p.m / self.g.D(x)**2
        return self.p.m * self.g.s(x)**2 / self.g.D(x)

    def substrate_phase(self, x):
        g, p, Q = self.g, self.p, self.Q
        if p.phase_mode == "proper_integral":
            return Q*g.ell(x), Q*g.s(x), Q*g.sp(x)
        if p.phase_mode == "literal_sr":
            phi = Q*(g.s(x)*x - float(g.s(p.R))*p.R)
            phip = Q*(g.s(x) + x*g.sp(x))
            phipp = Q*(2*g.sp(x) + x*(2*g.A/x**3))
            return phi, phip, phipp
        frozen_s = g.s(self.q)
        return Q*frozen_s*(x-p.R), Q*frozen_s, np.zeros_like(x)

    def potential_gradient_hessian(self, u):
        x = self.position(u)
        p, g, n = self.p, self.g, self.p.sites
        U = 0.0
        grad = np.zeros(n)
        H = np.zeros((n, n))
        if p.include_tidal:
            U += float(np.sum(p.m*p.c**2 * (g.D(x)-g.D(self.q)-g.Dp(self.q)*u)))
            grad += p.m*p.c**2 * (g.Dp(x)-g.Dp(self.q))
            H[np.arange(n), np.arange(n)] += p.m*p.c**2 * g.Dpp(x)

        phi, phip, phipp = self.substrate_phase(x)
        U += float(p.V0 * np.sum(g.D(x)*(1-np.cos(phi))))
        grad += p.V0 * (g.Dp(x)*(1-np.cos(phi)) + g.D(x)*np.sin(phi)*phip)
        H[np.arange(n), np.arange(n)] += p.V0 * (
            g.Dpp(x)*(1-np.cos(phi))
            + 2*g.Dp(x)*np.sin(phi)*phip
            + g.D(x)*(np.cos(phi)*phip**2 + np.sin(phi)*phipp))

        length = np.diff(g.ell(x))
        strain = length - self.l0
        mid = (x[1:] + x[:-1])/2
        d, dp, dpp = g.D(mid), g.Dp(mid), g.Dpp(mid)
        s, sp = g.s(x), g.sp(x)
        U += float(p.kappa/2 * np.sum(d * strain**2))
        for i, e in enumerate(strain):
            j = i + 1
            ei, ej = -s[i], s[j]
            eii, ejj = -sp[i], sp[j]
            db = d[i]
            dbi = dp[i]/2
            dbii = dpp[i]/4
            grad[i] += p.kappa/2 * (dbi*e**2 + 2*db*e*ei)
            grad[j] += p.kappa/2 * (dbi*e**2 + 2*db*e*ej)
            H[i,i] += p.kappa/2*(dbii*e**2 + 4*dbi*e*ei + 2*db*(ei**2 + e*eii))
            H[j,j] += p.kappa/2*(dbii*e**2 + 4*dbi*e*ej + 2*db*(ej**2 + e*ejj))
            off = p.kappa/2*(dbii*e**2 + 2*dbi*e*(ei+ej) + 2*db*ei*ej)
            H[i,j] += off
            H[j,i] += off
        return U, grad, H

    def potential(self, u):
        return self.potential_gradient_hessian(u)[0]

    def grad(self, u):
        return self.potential_gradient_hessian(u)[1]

    def hessian(self, u):
        return self.potential_gradient_hessian(u)[2]

    def equilibrium(self, starts=4):
        p = self.p
        best = None
        bound = 0.45 * p.dr
        bounds = [(-bound, bound)]*p.sites
        for seed in range(starts):
            rng = np.random.default_rng(seed)
            init = (np.zeros(p.sites) if seed == 0 else 0.08*p.dr*rng.standard_normal(p.sites))
            opt = minimize(self.potential, init, method="L-BFGS-B", jac=self.grad,
                           bounds=bounds, options={"maxiter": 3000, "gtol": 1e-9,
                                                   "ftol": 1e-13, "maxls": 50})
            # Deterministic Newton polishing of an interior, positive-Hessian
            # stationary solution. Never declare constrained convergence a
            # stationary equilibrium if a bound is active.
            x=opt.x.copy()
            for _ in range(8):
                uval, gu, Hu = self.potential_gradient_hessian(x)
                if np.max(np.abs(gu))<1e-11:
                    break
                if np.linalg.eigvalsh(Hu)[0] <= 0:
                    break
                step = np.linalg.solve(Hu, gu)
                accepted=False
                for alpha_ls in [1.0,0.5,0.25,0.125,0.0625]:
                    trial=x-alpha_ls*step
                    if np.max(np.abs(trial)) >= bound:
                        continue
                    if self.potential(trial)<uval:
                        x=trial
                        accepted=True
                        break
                if not accepted:
                    break
            gin=float(np.max(np.abs(self.grad(x))))
            cand = {"u":x, "energy":float(self.potential(x)),
                    "grad_inf":gin,"success":bool(opt.success),
                    "stationary":bool(gin<1e-6 and np.max(np.abs(x))<bound-1e-7)}
            if best is None or cand["energy"] < best["energy"]:
                best = cand
        return best

    def modes(self, u):
        H = self.hessian(u)
        mass = self.mass(u)
        w2, vec = eigh(H, np.diag(mass))
        return w2, vec, mass


def match_modes(vec_ref, vec_other, w2_ref, w2_other, count=20):
    """Mode assignment in a common radial-coordinate displacement basis.

    M-normalized modes from different metrics are rescaled to unit Euclidean
    length before a modal-assurance (MAC) / Hungarian assignment.
    """
    count = min(count, vec_ref.shape[1], vec_other.shape[1])
    # Keep all modes as assignment candidates; later report the lowest 20.
    af = vec_ref / np.maximum(np.linalg.norm(vec_ref, axis=0), 1e-300)
    bs = vec_other / np.maximum(np.linalg.norm(vec_other, axis=0), 1e-300)
    mac = np.abs(af.T @ bs)**2
    row, col = linear_sum_assignment(1-mac)
    inv = np.empty(vec_ref.shape[1], dtype=int)
    inv[row] = col
    ids = np.where(w2_ref > 1e-10)[0][:count]
    matched = inv[ids]
    overlap = mac[ids, matched]
    of = np.sqrt(np.maximum(w2_ref[ids], 0))
    os = np.sqrt(np.maximum(w2_other[matched], 0))
    alpha = float(np.dot(of,os)/np.dot(of,of))
    residual = float(np.max(np.abs(os-alpha*of)/np.maximum(alpha*of,1e-300)))
    trusted = overlap >= 0.8
    trusted_result = None
    if int(np.sum(trusted)) >= 5:
        aft = of[trusted]
        ast = os[trusted]
        a_trusted = float(np.dot(aft,ast)/np.dot(aft,aft))
        trusted_result = {"n":int(np.sum(trusted)), "alpha":a_trusted,
                          "residual":float(np.max(np.abs(ast-a_trusted*aft)/np.maximum(a_trusted*aft,1e-300)))}
    return {"alpha":alpha,"residual":residual,
            "high_mac_subset":trusted_result,"ref_modes":ids.tolist(),
            "mapped_modes":matched.tolist(),"mac":overlap.tolist(),
            "mac_min":float(np.min(overlap)),"mac_median":float(np.median(overlap)),
            "ambiguous_below_0_5":int(np.sum(overlap<0.5))}


def compare(model_ref, model_alt, count=20, starts=4):
    a=model_ref.equilibrium(starts)
    b=model_alt.equilibrium(starts)
    w2a,va,ma=model_ref.modes(a["u"])
    w2b,vb,mb=model_alt.modes(b["u"])
    sorted_a=np.sqrt(np.maximum(w2a,0))
    sorted_b=np.sqrt(np.maximum(w2b,0))
    good_a=sorted_a[sorted_a>1e-5][:count]
    good_b=sorted_b[sorted_b>1e-5][:count]
    if len(good_a)!=len(good_b):
        return {"physical_comparison_admissible":False,
                "reason":"insufficient positive modes in at least one branch",
                "positive_modes_ref":int(np.sum(w2a>1e-10)),
                "positive_modes_other":int(np.sum(w2b>1e-10)),
                "min_gen_eig_flat":float(w2a[0]),
                "min_gen_eig_other":float(w2b[0]),
                "flat_grad_inf":a["grad_inf"],
                "other_grad_inf":b["grad_inf"],
                "flat_stationary":a["stationary"],
                "other_stationary":b["stationary"]}
    alpha=float(np.dot(good_a,good_b)/np.dot(good_a,good_a))
    rank_res=float(np.max(np.abs(good_b-alpha*good_a)/(alpha*good_a)))
    m=match_modes(va,vb,w2a,w2b,count)
    return {"rank_alpha":alpha, "rank_residual":rank_res,
            "overlap_matching":m,
            "min_gen_eig_flat":float(w2a[0]),"min_gen_eig_other":float(w2b[0]),
            "min_mass_other":float(np.min(mb)),"max_mass_other":float(np.max(mb)),
            "flat_grad_inf":a["grad_inf"],"other_grad_inf":b["grad_inf"],
            "flat_success":a["success"],"other_success":b["success"],
            "flat_stationary":a["stationary"],"other_stationary":b["stationary"],
            "physical_comparison_admissible":bool(a["stationary"] and b["stationary"] and
                w2a[0]>-1e-8 and w2b[0]>-1e-8),
            "mode_identification_ambiguous":bool(m["mac_min"]<0.5)}


def run_probes(params=Parameters(), starts=4):
    base=SSZFKAction(params,flat=True)
    main=SSZFKAction(params,flat=False)
    result={"kind":"SSZ-FK proper-length and action consistency pilot",
            "scope":"EFFECTIVE_MODEL_ONLY; no astrophysical detection or QNM claims",
            "parameters":vars(params),"metric":{
                "D_R":float(main.g.D(params.R)),"s_R":float(main.g.s(params.R)),
                "proper_span":float(main.g.ell(params.R+(params.sites-1)*params.dr))},
            "action_ssz_vs_flat":compare(base,main,starts=starts)}
    probes = {
        "literal_sr_phase": {"phase_mode":"literal_sr"},
        "frozen_local_phase": {"phase_mode":"frozen_local"},
        "legacy_mass": {"mass_mode":"legacy_m_over_d2"},
        "constant_natural_proper_bonds": {"natural_bonds":"constant"},
        "no_substrate": {"V0":0.0},
        "no_substrate_no_tidal": {"V0":0.0,"include_tidal":False},
    }
    result["ablations"] = {}
    for name, overrides in probes.items():
        p2=Parameters(**{**vars(params),**overrides})
        result["ablations"][name]=compare(SSZFKAction(p2, flat=True),
                                          SSZFKAction(p2, flat=False),
                                          starts=(1 if p2.V0==0 else starts))
    return result


if __name__ == '__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=Path('ssz_fk_action_probe_v1.json'))
    p.add_argument('--sites',type=int,default=40)
    p.add_argument('--V0',type=float,default=0.0003)
    p.add_argument('--starts',type=int,default=4)
    a=p.parse_args()
    result=run_probes(Parameters(sites=a.sites,V0=a.V0), starts=a.starts)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in {'action_ssz_vs_flat','no_substrate_ssz_vs_flat'}},indent=2))
    print('action:', result['action_ssz_vs_flat']['rank_residual'],result['action_ssz_vs_flat']['overlap_matching']['residual'])
    for name, v in result['ablations'].items():
        print(name+':',round(v['rank_residual'],6) if 'rank_residual' in v else 'N/A',
              'MACmin',round(v['overlap_matching']['mac_min'],4) if 'overlap_matching' in v else 'N/A',
              'ADMISSIBLE' if v['physical_comparison_admissible'] else 'INVALID / OPEN')
