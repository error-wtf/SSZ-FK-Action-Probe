# ChatGPT Round II — Ordering-Constraint Probe (2026-10-11)

**Auftrag (ChatGPT, Auswertung von Bericht II):** Ersetze die künstliche
Einzelschranke `|u_i| ≤ 0.45·dr` durch echte Nachbar-Ordnungsbedingungen
`r_{i+1} − r_i ≥ ε`, reproduziere die fünf problematischen Intervalle und
prüfe erst dann, ob der behauptete Fold existiert.

**Korrektheit der Kritik:** Vollauf bestätigt. Die Box ist hinreichend, aber
nicht notwendig für Ordnung (`u_{i+1} − u_i > −Δr` ist die echte Bedingung);
ein KKT-Randhit allein beweist keinen Fold (Gegenbeispiel `U=(u−μ)²/2,
|u|≤a`).

## Methode

Werkzeug: `tools/ordering_constraint_probe_v1.py` (L-BFGS-B in
Gap-Variablen `v = (u_0, g_1, …, g_{n-1})`, `u_i = u_0 + Σ_{j≤i} g_j` — die
Ordnungsbedingung wird zur simplen Box in v; Gradient/Hessian via
Kettenregel, finite-difference-verifiziert mit max. rel. Fehler
2.1e-8 (Grad) / 3.6e-7 (Hessian)). Multistart (warm + 4 random),
KKT-Klassifikation; dazu **unbeschränkte Newton-Kontinuation** (keine
Constraints, nur Backtracking-Damping) mit λ_min-Monitoring und Bisektion
der Nulldurchgänge. `ε ∈ {0.05, 0.2, 0.5}·dr`, Sites N=40, seed 1104.

Artefakt: `SSZ_FK_ORDERING_CONSTRAINT_V1.json`
(sha256 `e3291769…`, strict JSON, Exit 0).
Figuren: `figures/fig3_ordering_constraint.png`,
`figures/anim_unconstrained_continuation.gif` (echte Live-Läufe).

## Ergebnis 1 — Die Box-Lücken sind Ordnungs-Artefakte

**Alle fünf problematischen Intervalle sind unter der echten
Ordnungsbedingung vollständig admissibel — bei allen drei ε:**

| Intervall | Box-Setup | ε=0.05·dr | ε=0.2·dr | ε=0.5·dr |
|---|---|---|---|---|
| gap1 0.2795–0.8998 | inadmissibel | **15/15 admissibel** | 12/12+3=15/15 | 15/15 |
| gap2 1.8580–2.3498 | inadmissibel | **12/12 admissibel** | 12/12 | 12/12 |
| gap3 2.6057–2.9502 | inadmissibel | **9/9 admissibel** | 9/9 | 9/9 |
| endzone 3.35–3.45 | fragmentiert | **5/5 admissibel** | 5/5 | 5/5 |
| endzone 4.05–4.10 | fragmentiert | **4/4 admissibel** | 4/4 | 4/4 |

Keine Bisektionskanten nötig (0 Admissibilitätswechsel in 45 Punkten × 3 ε).
Die Lösungen sind **streng innere Punkte**: minimale radiale Lücke
0.0409–0.0471 über alle unbeschränkten Lösungen (dr = 0.05), d.h. die
Ordnung `r_{i+1} > r_i` ist nirgends aktiv oder auch nur nahe aktiv. Es gab
niemals ein physikalisches Ordnungsproblem in diesen Bereichen — nur die
zu enge Verschiebungsbox.

**Konsequenz:** Die CONSTRAINT_BOUNDARY-Klassifikation von Bericht II in den
drei Lücken war ein **Modellierungsartefakt der Box**, kein physikalisches
Ereignis. ChatGPTs Einwand trifft exakt zu. Bericht II wird hiermit
korrigiert: "kein inneres Gleichgewicht in den Lücken" ist **falsch** — es
gibt eins, die Box hat es nur verboten.

## Ergebnis 2 — Der unbeschränkte λ_min-Verlauf (nuanciert, ehrlich)

Unbeschränkte Newton-Kontinuation (warm gestartet, 89 Punkte, Δr_s=0.05):
43/89 konvergiert. λ_min > 0 für rs 0.00–0.95, **Vorzeichenwechsel bei
rs ≈ 0.97** (8.1e-3 bei 0.95 → −2.3e-3 bei 1.00), danach Stall-Bereich bis
~1.9, wieder konvergiert 1.9–3.05, wieder negativ ab ~3.10.

**Aber:** Die Stall-Region ist **pfadabhängig** — beim Direkteinstieg bei
rs=0.90 (warm von 0.0 mit grober Schrittweite) stagniert Newton bei
grad_inf ≈ 7e-2, λ_min ≈ −1.4; derselbe rs-Punkt konvergierte im
0.05er-Lauf. Das ist konsistent mit mehreren unbeschränkten Branches in
dieser Region (mögliche Sliding-Mode-Struktur, vgl. Sokoloff et al. 1978,
unpinned phase), ist aber **kein sauber lokalisierter Fold**: dafür wäre
Pseudo-Arclength-Kontinuation mit Branch-Tracking nötig. Der λ_min-Wechsel
bei rs ≈ 0.97 ist ein **Fold-Kandidat**, kein belegter Fold.

## Status der ChatGPT-Prüfliste

| Punkt | Status |
|---|---|
| KKT-Randlösung ≠ physikalischer Fold | **bestätigt** — Box-Artefakt nachgewiesen |
| Ordnungsbedingung statt Box | **erledigt** — alle 5 Intervalle admissibel |
| Fold-Behauptung | **zurückgezogen** für die 3 Lücken; Fold-Kandidat rs≈0.97 offen |
| Energie/Constraints/Gradient/Hessian/Modenresiduen an denselben Punkten | Energie + Gradient + Hessian erledigt; Modenresiduen offen |
| Kontinuation der unbeschränkten Gleichungen + Verschwindende Eigenwerte + Branch-Umkehrpunkte | Vorläufig erledigt (λ≈0 bei rs≈0.97); Branch-Tracking (Pseudo-Arclength) offen |
| Fixed-Domain-Konvergenz | offen (Runde I, Punkt 3) |
| SSZ-vs-GR-Nullvergleich | offen (Runde I, Punkt 4) |

## Nächste Schritte (Priorität)

1. **Pseudo-Arclength-Kontinuation** der unbeschränkten stationären
   Gleichungen mit explizitem Branch-Tracking durch rs ∈ [0.9, 1.2] —
   Fold-Kandidat rs ≈ 0.97 verifizieren oder verwerfen (Hysterese-Test
   aufwärts/abwärts inklusive).
2. Modenresiduen der neuen (unbeschränkten) Lösungen an denselben
   rs-Punkten — Anschluss an die MAC-Diagnostik aus fig2.
3. Erst danach: Fixed-Domain-Konvergenz und SSZ-vs-GR-Nullvergleich
   (ChatGPT: "Wenn die Struktur diesen Test übersteht …").

## Reproduktion

```bash
PY=/home/error/physics/venv/bin/python
$PY tools/ordering_constraint_probe_v1.py \
    --output SSZ_FK_ORDERING_CONSTRAINT_V1.json   # ~90 s
$PY tools/make_ordering_figures_v1.py              # ~60 s (GIF live-Läufe)
```

Umgebung wie im Full-Pipeline-Report (Kali, Python 3.14.7 venv,
OMP_NUM_THREADS=2).
