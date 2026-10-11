import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from ssz_fk_action_probe_v1 import Metric, Parameters, SSZFKAction, match_modes


def model(flat=False, **kw):
    return SSZFKAction(replace(Parameters(sites=12), **kw), flat=flat)


def test_metric_proper_length_and_phase_derivative():
    m=Metric()
    x=np.linspace(12,14,21)
    assert np.max(np.abs(m.s(x)*m.D(x)-1))<1e-14
    assert np.max(np.abs(np.gradient(m.ell(x),x,edge_order=2)[1:-1]-m.s(x[1:-1])))<2e-5
    mdl=model()
    for mode in ('proper_integral','literal_sr','frozen_local'):
        mdl=model(phase_mode=mode)
        phase, deriv, second=mdl.substrate_phase(mdl.q)
        h=1e-5
        plus=mdl.substrate_phase(mdl.q+h)[0]
        minus=mdl.substrate_phase(mdl.q-h)[0]
        assert np.max(np.abs((plus-minus)/(2*h)-deriv))<1e-7
        dplus=mdl.substrate_phase(mdl.q+h)[1]
        dminus=mdl.substrate_phase(mdl.q-h)[1]
        assert np.max(np.abs((dplus-dminus)/(2*h)-second))<1e-7
    assert np.max(np.abs(model(phase_mode='literal_sr').substrate_phase(mdl.q)[1]-mdl.Q))<1e-12


@pytest.mark.parametrize('flat', [False, True])
@pytest.mark.parametrize('phase', ['proper_integral','literal_sr','frozen_local'])
def test_analytic_grad_hessian(flat,phase):
    mdl=model(flat=flat,phase_mode=phase)
    rng=np.random.default_rng(20261011)
    u=0.01*mdl.p.dr*rng.normal(size=mdl.p.sites)
    en,g,hess=mdl.potential_gradient_hessian(u)
    delta=1e-6
    gn=np.empty(mdl.p.sites)
    Hn=np.empty_like(hess)
    for i in range(mdl.p.sites):
        step=np.zeros(mdl.p.sites);step[i]=delta
        gn[i]=(mdl.potential(u+step)-mdl.potential(u-step))/(2*delta)
        Hn[:,i]=(mdl.grad(u+step)-mdl.grad(u-step))/(2*delta)
    assert np.max(np.abs(g-gn))<3e-8
    assert np.max(np.abs(hess-Hn))<1e-5
    assert np.max(np.abs(hess-hess.T))<1e-12


def test_mass_from_point_particle_action():
    mdl=model()
    x=mdl.q+0.001*mdl.p.dr
    D=mdl.g.D(x);s=mdl.g.s(x)
    c=mdl.p.c;m=mdl.p.m
    v=1e-3*c*D/s
    E_rel=-m*c*c*D*np.sqrt(1-s*s*v*v/(c*c*D*D))
    E_rest=-m*c*c*D
    quadratic=0.5*m*s*s/D*v*v
    assert np.max(np.abs((E_rel-E_rest-quadratic)/quadratic))<1e-6
    assert np.max(np.abs(mdl.mass(np.zeros(mdl.p.sites))-m*s*s/D))<0.01


def test_flat_zero_rs_and_free_chain():
    mdl=model(flat=True,V0=0.,include_tidal=False)
    z=np.zeros(mdl.p.sites)
    _,grad,H=mdl.potential_gradient_hessian(z)
    mass=mdl.mass(z)
    assert np.allclose(mass,1.0)
    assert np.max(np.abs(grad))<1e-10
    vals, vec, _=mdl.modes(z)
    ana=2*np.sqrt(mdl.p.kappa/mdl.p.m)*np.sin(np.arange(mdl.p.sites)*np.pi/(2*mdl.p.sites))
    assert np.max(np.abs(np.sqrt(np.maximum(vals,0))[1:]-ana[1:]))<1e-11
    assert np.max(np.abs(vec.T@np.diag(mass)@vec-np.eye(mdl.p.sites)))<1e-11


def test_mass_normalized_resolvent_and_sum_rule():
    mdl=model(V0=0.0003)
    ue=mdl.equilibrium(starts=2)['u']
    vals, vec, mass=mdl.modes(ue)
    H=mdl.hessian(ue)
    b=np.zeros(mdl.p.sites);b[2]=1
    d=np.zeros(mdl.p.sites);d[7]=1
    z=(3+0.3j)**2
    direct=d@np.linalg.solve(H-z*np.diag(mass),b)
    spec=np.sum((vec.T@d)*(vec.T@b)/(vals-z))
    assert abs(spec-direct)<1e-10
    assert abs(np.sum((vec.T@b)**2)-b@np.diag(1/mass)@b)<1e-12


def test_assignment_is_not_rank_matching():
    v=np.eye(4)
    swapped=v[:,[1,0,3,2]]
    a=np.array([1.,4.,9.,16.])
    b=a[[1,0,3,2]]
    m=match_modes(v,swapped,a,b,count=4)
    assert m['mapped_modes']==[1,0,3,2]
    assert m['residual']<1e-13
    assert m['mac_min']>0.999999


def test_exact_flat_branch_agrees_with_zero_rs():
    p=Parameters(sites=12,r_s=0)
    a=SSZFKAction(p,flat=False)
    b=SSZFKAction(p,flat=True)
    u=np.linspace(-1e-3,1e-3,p.sites)
    for x,y in zip(a.potential_gradient_hessian(u),b.potential_gradient_hessian(u)):
        assert np.allclose(x,y,atol=1e-12)
    assert np.allclose(a.mass(u),b.mass(u))
