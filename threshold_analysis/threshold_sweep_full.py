"""
Threshold sensitivity sweep, Case 1 (double-knockout, iAF1260, anaerobic succinate).

Re-uses the exact RapidKnock pipeline (P1/P2 FVA candidate screening, lethality
pruning, combinatorial testing) at 5 flux-threshold values, keeping every other
parameter identical to the reported Case 1 run. Records, per threshold:
  A) candidate pool size (reactions passing the P1/P2 FVA screen)
  B) wall-clock runtime
  C) recovered succinate yield of the best viable combination
"""

import cobra
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
import numpy as np
from itertools import combinations
import warnings
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

# ── Fixed Case 1 conditions (unchanged from RapidKnock.py) ─────────────────
MODEL_FILE        = "iAF1260.xml"
BIOMASS_RXN       = "BIOMASS_Ec_iAF1260_core_59p81M"
SUCCINATE_EX      = "EX_succ_e"
GLUCOSE_EX        = "EX_glc__D_e"
OXYGEN_EX         = "EX_o2_e"
SOLVER            = "glpk"
ANAEROBIC         = True
MIN_BIOMASS_FRAC  = 0.01
MAX_KO            = 2
GREEDY_STEP_SIZE  = 2
MIN_SUCC_YIELD    = 0.50
N_PROCESSES       = -1
TURN_OFF_REACTIONS = []
EXCLUDED_REACTIONS = ["ATPM"]

THRESHOLDS = [1e-5, 1e-6, 1e-7, 1e-8]

# ── Helpers (identical to RapidKnock.py) ────────────────────────────────────

def is_transport(rxn):
    base_ids = set(m.id.rsplit("_", 1)[0] for m in rxn.metabolites.keys())
    return len(base_ids) == 1

def get_rxns_with_gpr(m, locked_kos):
    return set(
        r.id for r in m.reactions
        if r.gene_reaction_rule.strip() != ""
        and not (r.lower_bound == 0 and r.upper_bound == 0)
        and not is_transport(r)
        and r.id not in EXCLUDED_REACTIONS
        and r.id not in TURN_OFF_REACTIONS
        and r.id not in locked_kos
    )

def run_fva_pipeline(m, flux_threshold_frac):
    with m:
        m.objective = BIOMASS_RXN
        sol = m.optimize()
        ctx_biomass = sol.objective_value if sol.status == "optimal" else 0.0
        fva_p1 = flux_variability_analysis(m, fraction_of_optimum=1.0, processes=N_PROCESSES)

    with m:
        m.reactions.get_by_id(BIOMASS_RXN).lower_bound = MIN_BIOMASS_FRAC * ctx_biomass
        m.objective = SUCCINATE_EX
        p2_sol = m.optimize()
        fva_p2 = flux_variability_analysis(m, fraction_of_optimum=1.0, processes=N_PROCESSES)

    glc_uptake     = abs(fva_p1.loc[GLUCOSE_EX, "minimum"])
    flux_threshold = flux_threshold_frac * glc_uptake

    locked_kos = [r.id for r in m.reactions
                  if r.lower_bound == 0 and r.upper_bound == 0
                  and r.id not in EXCLUDED_REACTIONS]
    rxns_with_gpr = get_rxns_with_gpr(m, locked_kos)

    results = []
    for rxn_id in fva_p1.index:
        p1_min = fva_p1.loc[rxn_id, "minimum"]
        p2_min = fva_p2.loc[rxn_id, "minimum"]
        p2_max = fva_p2.loc[rxn_id, "maximum"]

        is_candidate = (
            rxn_id in rxns_with_gpr
            and abs(p1_min) > flux_threshold
            and (p2_max < p1_min or abs(p2_min) < abs(p1_min))
        )
        results.append({"reaction": rxn_id, "p1_min": p1_min, "is_candidate": is_candidate})

    df = pd.DataFrame(results)
    knockout_candidates = df[
        df["is_candidate"]
        & ~df["reaction"].isin(TURN_OFF_REACTIONS)
        & ~df["reaction"].isin(EXCLUDED_REACTIONS)
    ].sort_values("p1_min", ascending=False)

    return knockout_candidates

def simulate_kos(base_model, ko_set):
    with base_model:
        for rxn_id in ko_set:
            base_model.reactions.get_by_id(rxn_id).knock_out()
        base_model.objective = BIOMASS_RXN
        sol = base_model.optimize()
        if sol.status != "optimal":
            return 0.0, 0.0, 0.0, sol.status
        growth    = sol.objective_value
        succ_flux = sol.fluxes.get(SUCCINATE_EX, 0.0)
        glc_flux  = abs(sol.fluxes.get(GLUCOSE_EX))
        return growth, succ_flux, round(succ_flux / glc_flux, 6), sol.status

def run_case1_at_threshold(flux_threshold_frac):
    t0 = time.time()

    model = cobra.io.read_sbml_model(MODEL_FILE)
    model.solver = SOLVER
    if ANAEROBIC:
        model.reactions.get_by_id(OXYGEN_EX).lower_bound = 0
    for rxn_id in TURN_OFF_REACTIONS:
        model.reactions.get_by_id(rxn_id).knock_out()

    with model:
        model.objective = BIOMASS_RXN
        wt_biomass = model.optimize().objective_value

    working_model = model.copy()
    ko_candidates_df = run_fva_pipeline(working_model, flux_threshold_frac)
    candidate_pool_size = len(ko_candidates_df)

    if ko_candidates_df.empty:
        return {"threshold": flux_threshold_frac, "candidate_pool_size": 0,
                "runtime_s": round(time.time() - t0, 2), "best_yield": None}

    candidates = ko_candidates_df["reaction"].tolist()
    lethal_rxns = set()
    for rxn_id in candidates:
        growth, _, _, _ = simulate_kos(working_model, [rxn_id])
        if growth <= MIN_BIOMASS_FRAC * wt_biomass:
            lethal_rxns.add(rxn_id)
    viable_candidates = [r for r in candidates if r not in lethal_rxns]

    if len(viable_candidates) < GREEDY_STEP_SIZE:
        return {"threshold": flux_threshold_frac, "candidate_pool_size": candidate_pool_size,
                "runtime_s": round(time.time() - t0, 2), "best_yield": None}

    combos = list(combinations(viable_candidates, GREEDY_STEP_SIZE))
    best_yield = None
    for combo in combos:
        growth, succ_flux, succ_yield, status = simulate_kos(working_model, list(combo))
        if growth > MIN_BIOMASS_FRAC * wt_biomass:
            if best_yield is None or succ_yield > best_yield:
                best_yield = succ_yield

    runtime_s = round(time.time() - t0, 2)
    return {"threshold": flux_threshold_frac, "candidate_pool_size": candidate_pool_size,
            "runtime_s": runtime_s, "best_yield": best_yield}

# ── Run sweep ─────────────────────────────────────────────────────────────

rows = []
for thr in THRESHOLDS:
    print(f"Running threshold = {thr:.0e} ...")
    row = run_case1_at_threshold(thr)
    rows.append(row)
    print(f"  candidates={row['candidate_pool_size']}, "
          f"runtime={row['runtime_s']}s, best_yield={row['best_yield']}")

df = pd.DataFrame(rows)
df.to_excel("threshold_sweep_case1.xlsx", index=False)
print("\nSaved threshold_sweep_case1.xlsx")
print(df.to_string(index=False))

# ── Plot 3-panel figure directly from the sweep results ────────────────────

df = df.sort_values("threshold")
x = df["threshold"].values
labels = [f"1e-{int(-np.log10(t))}" for t in x]
c_main = '#1b6ca8'

def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.4)
    ax.spines['bottom'].set_linewidth(1.4)
    ax.tick_params(axis='both', width=1.4, labelsize=11)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight('bold')

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
panels = [
    ("A", df["candidate_pool_size"], "Candidate pool size (# reactions)"),
    ("B", df["runtime_s"], "Runtime (s)"),
    ("C", df["best_yield"], "Recovered yield (mol succinate / mol glucose)"),
]

for ax, (label, y, ylab) in zip(axes, panels):
    ax.plot(range(len(x)), y, 'o-', color=c_main, linewidth=2.5, markersize=9)
    ax.set_xticks(range(len(x)))
    ax.set_xticklabels(labels, fontsize=11, fontweight='bold')
    ax.set_xlabel('Flux threshold (fraction of carbon uptake)', fontsize=12, fontweight='bold')
    ax.set_ylabel(ylab, fontsize=11, fontweight='bold')
    ax.set_title(label, loc='left', fontsize=16, fontweight='bold')
    style_ax(ax)

plt.tight_layout()
plt.savefig('threshold_sweep_case1.png', dpi=300, bbox_inches='tight')
print("Saved threshold_sweep_case1.png")
