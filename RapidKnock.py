"""
FVA-based Greedy Iterative KO Pipeline for Succinate Overproduction

Each round:
  1. Run P1 + P2 FVA on current KO background
  2. Derive knockout candidates from P1/P2 flux comparison
  3. Prune lethals within current KO background
  4. Test all combos of GREEDY_STEP_SIZE, pick best viable by succ_yield
  5. Lock chosen KOs, repeat until MAX_KO reached
Then:
  - Collect all simulations across all rounds
  - Report top results meeting yield and biomass criteria
  - Save all results and timing to Excel
"""

import cobra
from cobra.flux_analysis import flux_variability_analysis
import pandas as pd
from itertools import combinations
import warnings
import time
warnings.filterwarnings("ignore")

# ── Parameters ─────────────────────────────────────────────────────────────────

MODEL_FILE        = "iAF1260.xml"
BIOMASS_RXN       = "BIOMASS_Ec_iAF1260_core_59p81M"
SUCCINATE_EX      = "EX_succ_e"
GLUCOSE_EX        = "EX_glc__D_e"
OXYGEN_EX         = "EX_o2_e"
SOLVER            = "glpk"
ANAEROBIC         = True
MIN_BIOMASS_FRAC  = 0.01
FLUX_THRESHOLD    = 0.000000001
MAX_KO            = 4
GREEDY_STEP_SIZE  = 4
MIN_SUCC_YIELD    = 0.50
N_PROCESSES       = -1
TURN_OFF_REACTIONS = []
EXCLUDED_REACTIONS = ["ATPM"]

# ──────────────────────────────────────────────────────────────────────────────

t_start = time.time()

# ── Helpers ────────────────────────────────────────────────────────────────────

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

def run_fva_pipeline(m, wt_biomass):
    with m:
        m.objective = BIOMASS_RXN
        sol = m.optimize()
        ctx_biomass = sol.objective_value if sol.status == "optimal" else 0.0
        fva_p1 = flux_variability_analysis(m, fraction_of_optimum=1.0, processes=N_PROCESSES)

    with m:
        m.reactions.get_by_id(BIOMASS_RXN).lower_bound = MIN_BIOMASS_FRAC * ctx_biomass
        m.objective = SUCCINATE_EX
        p2_sol = m.optimize()
        p2_succ = p2_sol.objective_value if p2_sol.status == "optimal" else 0.0
        fva_p2 = flux_variability_analysis(m, fraction_of_optimum=1.0, processes=N_PROCESSES)

    glc_uptake     = abs(fva_p1.loc[GLUCOSE_EX, "minimum"])
    flux_threshold = FLUX_THRESHOLD * glc_uptake

    locked_kos = [r.id for r in m.reactions
                  if r.lower_bound == 0 and r.upper_bound == 0
                  and r.id not in EXCLUDED_REACTIONS]
    rxns_with_gpr = get_rxns_with_gpr(m, locked_kos)

    results = []
    for rxn_id in fva_p1.index:
        p1_min = fva_p1.loc[rxn_id, "minimum"]
        p1_max = fva_p1.loc[rxn_id, "maximum"]
        p2_min = fva_p2.loc[rxn_id, "minimum"]
        p2_max = fva_p2.loc[rxn_id, "maximum"]

        is_candidate = (
            rxn_id in rxns_with_gpr
            and abs(p1_min) > flux_threshold
            and (
                p2_max < p1_min
                or abs(p2_min) < abs(p1_min)
            )
        )
        results.append({
            "reaction"     : rxn_id,
            "p1_min"       : round(p1_min, 4),
            "p1_max"       : round(p1_max, 4),
            "p2_min"       : round(p2_min, 4),
            "p2_max"       : round(p2_max, 4),
            "is_candidate" : is_candidate,
        })

    df = pd.DataFrame(results)
    knockout_candidates = df[
        df["is_candidate"]
        & ~df["reaction"].isin(TURN_OFF_REACTIONS)
        & ~df["reaction"].isin(EXCLUDED_REACTIONS)
    ].sort_values("p1_min", ascending=False)

    return knockout_candidates, ctx_biomass, glc_uptake, p2_succ

def simulate_kos(base_model, ko_set, wt_biomass):
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

# ── 1. Load model ──────────────────────────────────────────────────────────────

model = cobra.io.read_sbml_model(MODEL_FILE)
model.solver = SOLVER

if ANAEROBIC:
    model.reactions.get_by_id(OXYGEN_EX).lower_bound = 0

for rxn_id in TURN_OFF_REACTIONS:
    model.reactions.get_by_id(rxn_id).knock_out()
    print(f"  Turned off: {rxn_id}")

print("Model loaded:")
print(f"  Reactions  : {len(model.reactions)}")
print(f"  Metabolites: {len(model.metabolites)}")
print(f"  Genes      : {len(model.genes)}")
print(f"  Solver     : {SOLVER}")
print(f"  Anaerobic  : {ANAEROBIC}")
print(f"  MAX_KO     : {MAX_KO}  (step size: {GREEDY_STEP_SIZE})")
print(f"  Processes  : {N_PROCESSES}")

# ── 2. WT biomass ──────────────────────────────────────────────────────────────

with model:
    model.objective = BIOMASS_RXN
    wt_biomass = model.optimize().objective_value
print(f"\n  Base WT max biomass: {wt_biomass:.4f}")

# ── 3. Greedy Iterative KO Loop ────────────────────────────────────────────────

rounds_schedule = []
remaining_budget = MAX_KO
while remaining_budget > 0:
    rounds_schedule.append(min(GREEDY_STEP_SIZE, remaining_budget))
    remaining_budget -= rounds_schedule[-1]

print(f"\n--- Greedy Iterative KO Selection (MAX_KO={MAX_KO}, step={GREEDY_STEP_SIZE}) ---")
print(f"  Round schedule: {rounds_schedule}")

selected_kos   = []
round_records  = []
all_round_sims = []

for round_idx, step in enumerate(rounds_schedule, start=1):
    print(f"\n{'='*60}")
    print(f"  ROUND {round_idx}/{len(rounds_schedule)} - selecting {step} KO(s)")
    print(f"  Locked KOs so far: {selected_kos or 'none'}")

    working_model = model.copy()
    for rxn_id in selected_kos:
        working_model.reactions.get_by_id(rxn_id).knock_out()

    print(f"  Running P1+P2 FVA...")
    ko_candidates_df, ctx_biomass, glc_uptake, p2_succ = run_fva_pipeline(
        working_model, wt_biomass
    )

    print(f"    Biomass in context : {ctx_biomass:.4f}")
    print(f"    Glucose uptake     : {glc_uptake:.4f}")
    print(f"    Max succinate (P2) : {p2_succ:.4f}")
    print(f"    Candidates         : {len(ko_candidates_df)}")

    if not ko_candidates_df.empty:
        print(ko_candidates_df[["reaction","p1_min","p1_max",
                                "p2_min","p2_max"]].to_string(index=False))
        print("\n    Reaction details:")
        for rxn_id in ko_candidates_df["reaction"].tolist():
            rxn = working_model.reactions.get_by_id(rxn_id)
            print(f"      {rxn_id:20s}: {rxn.reaction[:55]:55s} "
                  f"GPR={rxn.gene_reaction_rule[:30]}")

    if ko_candidates_df.empty:
        print("    No candidates found - stopping greedy search.")
        break

    candidates  = ko_candidates_df["reaction"].tolist()
    lethal_rxns = set()
    for rxn_id in candidates:
        growth, _, _, _ = simulate_kos(working_model, [rxn_id], wt_biomass)
        if growth <= MIN_BIOMASS_FRAC * wt_biomass:
            lethal_rxns.add(rxn_id)

    viable_candidates = [r for r in candidates if r not in lethal_rxns]
    print(f"\n    Pruned {len(lethal_rxns)} lethals")
    print(f"    Viable candidates : {len(viable_candidates)}")

    if len(viable_candidates) < step:
        print(f"    Not enough viable candidates ({len(viable_candidates)}) "
              f"for step size {step}. Stopping.")
        break

    combos = list(combinations(viable_candidates, step))
    print(f"    Testing {len(combos)} combination(s) of size {step}...")

    best_record = None

    for combo in combos:
        growth, succ_flux, succ_yield, status = simulate_kos(
            working_model, list(combo), wt_biomass
        )
        viable   = growth > MIN_BIOMASS_FRAC * wt_biomass
        full_kos = selected_kos + list(combo)

        record = {
            "round"          : round_idx,
            "round_combo"    : ", ".join(combo),
            "cumulative_kos" : ", ".join(full_kos),
            "n_kos"          : len(full_kos),
            "growth"         : round(growth, 4),
            "succ_flux"      : round(succ_flux, 4),
            "succ_yield"     : round(succ_yield, 6),
            "viable"         : viable,
        }
        all_round_sims.append(record)

        if viable and (best_record is None or succ_yield > best_record["succ_yield"]):
            best_record = record

    if best_record is None:
        print(f"    No viable combination in round {round_idx}. Stopping.")
        break

    chosen = best_record["round_combo"].split(", ")
    selected_kos.extend(chosen)
    round_records.append(best_record)

    print(f"\n    Best combo  : {best_record['round_combo']}")
    print(f"    Cumulative  : {best_record['cumulative_kos']}")
    print(f"    Growth      : {best_record['growth']:.4f}  "
          f"Succ flux: {best_record['succ_flux']:.4f}  "
          f"Yield: {best_record['succ_yield']:.4f}")

# ── 4. Round-by-round summary ──────────────────────────────────────────────────

print("\n--- Round-by-Round Summary ---")
if round_records:
    rdf = pd.DataFrame(round_records)
    print(rdf[["round","round_combo","cumulative_kos","n_kos",
               "growth","succ_flux","succ_yield"]].to_string(index=False))

# ── 5. Final result ────────────────────────────────────────────────────────────

print(f"\n{'='*60}")
print("--- Final Result ---")

total_time = round(time.time() - t_start, 2)

sim_df  = pd.DataFrame(all_round_sims)
passing = sim_df[
    sim_df["viable"]
    & (sim_df["succ_yield"] >= MIN_SUCC_YIELD)
].sort_values("succ_yield", ascending=False)

if not passing.empty:
    print(f"  {len(passing)} result(s) met all criteria across all rounds.")
    print(f"\n  Top results:")
    print(passing[["cumulative_kos","n_kos","growth",
                   "succ_flux","succ_yield"]].head(10).to_string(index=False))

    best_row        = passing.iloc[0]
    best_growth     = best_row["growth"]
    best_succ_yield = best_row["succ_yield"]
    best_succ_flux  = best_row["succ_flux"]
    reported_kos    = best_row["cumulative_kos"].split(", ")

    print(f"\n  Best result:")
    print(f"  Reported KO set : {reported_kos}")
    print(f"  Growth          : {best_growth:.4f}  "
          f"(min required: {MIN_BIOMASS_FRAC * wt_biomass:.4f})")
    print(f"  Succinate flux  : {best_succ_flux:.4f}")
    print(f"  Succinate yield : {best_succ_yield:.4f}  "
          f"(min required: {MIN_SUCC_YIELD})")

    all_results = passing.copy()
    if TURN_OFF_REACTIONS:
        prefix = ", ".join(TURN_OFF_REACTIONS) + ", "
        all_results["all_kos"]      = prefix + all_results["cumulative_kos"]
        all_results["n_kos_total"]  = all_results["n_kos"] + len(TURN_OFF_REACTIONS)
    else:
        all_results["all_kos"]      = all_results["cumulative_kos"]
        all_results["n_kos_total"]  = all_results["n_kos"]

    export_df = all_results[["all_kos","n_kos_total","growth",
                              "succ_flux","succ_yield"]].copy()
    export_df.columns = ["all_knockouts","n_kos","growth","succ_flux","succ_yield"]

else:
    print(f"  No result met criteria: viable, succ_yield >= {MIN_SUCC_YIELD}")
    print(f"  Reporting best viable result across all rounds.")
    viable_all = sim_df[sim_df["viable"]].sort_values("succ_yield", ascending=False)
    if not viable_all.empty:
        best_row        = viable_all.iloc[0]
        best_growth     = best_row["growth"]
        best_succ_yield = best_row["succ_yield"]
        best_succ_flux  = best_row["succ_flux"]
        reported_kos    = best_row["cumulative_kos"].split(", ")

        all_results = viable_all.copy()
        if TURN_OFF_REACTIONS:
            prefix = ", ".join(TURN_OFF_REACTIONS) + ", "
            all_results["all_kos"]     = prefix + all_results["cumulative_kos"]
            all_results["n_kos_total"] = all_results["n_kos"] + len(TURN_OFF_REACTIONS)
        else:
            all_results["all_kos"]     = all_results["cumulative_kos"]
            all_results["n_kos_total"] = all_results["n_kos"]

        export_df = all_results[["all_kos","n_kos_total","growth",
                                  "succ_flux","succ_yield"]].copy()
        export_df.columns = ["all_knockouts","n_kos","growth","succ_flux","succ_yield"]

        print(f"\n  Reported KO set : {reported_kos}")
        print(f"  Growth          : {best_growth:.4f}")
        print(f"  Succinate flux  : {best_succ_flux:.4f}")
        print(f"  Succinate yield : {best_succ_yield:.4f}")
    else:
        print("  No viable result found at all.")
        best_growth = best_succ_yield = best_succ_flux = 0.0
        reported_kos = []
        export_df = pd.DataFrame()

# ── 6. Save to Excel ───────────────────────────────────────────────────────────

out_file = "RapidKnock_results.xlsx"
with pd.ExcelWriter(out_file, engine="openpyxl") as writer:

    if not export_df.empty:
        export_df.to_excel(writer, sheet_name="all_results", index=False)

    summary_df = pd.DataFrame({
        "parameter" : [
            "model_file",
            "solver",
            "anaerobic",
            "min_biomass_frac",
            "flux_threshold",
            "max_ko",
            "greedy_step_size",
            "min_succ_yield",
            "turn_off_reactions",
            "wt_biomass",
            "total_time_s",
        ],
        "value" : [
            MODEL_FILE,
            SOLVER,
            ANAEROBIC,
            MIN_BIOMASS_FRAC,
            FLUX_THRESHOLD,
            MAX_KO,
            GREEDY_STEP_SIZE,
            MIN_SUCC_YIELD,
            ", ".join(TURN_OFF_REACTIONS) if TURN_OFF_REACTIONS else "none",
            round(wt_biomass, 4),
            total_time,
        ]
    })
    summary_df.to_excel(writer, sheet_name="run_summary", index=False)

print(f"\n  Results saved to: {out_file}")
print(f"\n  Total time: {total_time}s")