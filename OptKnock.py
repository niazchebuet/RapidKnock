import cobra
import straindesign as sd
import pandas as pd
import time

def is_knockout_value(val):
    if isinstance(val, str):
        return val.strip().upper() in {"KO", "KNOCKOUT", "DEL", "DELETE"}
    try:
        return val < 0
    except Exception:
        return False

def main():
    t_start = time.time()

    # ----------------------
    # Load model
    # ----------------------
    model = cobra.io.read_sbml_model("iAF1260.xml")

    # ----------------------
    # Anaerobic condition
    # ----------------------
    model.reactions.get_by_id("EX_o2_e").lower_bound = 0

    # ----------------------
    # Pre-fixed knockouts — must match TURN_OFF_REACTIONS in GreedyKnock
    # ----------------------
    turn_off_reactions = []
    for rxn_id in turn_off_reactions:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).knock_out()
            print(f"  Turned off: {rxn_id}")

    # ----------------------
    # EXCLUSION LIST
    # ----------------------
    excluded_rxns = ["ATPM"]

    # ----------------------
    # Biomass & target
    # ----------------------
    biomass_rxn = "BIOMASS_Ec_iAF1260_core_59p81M"
    target_rxn  = "EX_succ_e"
    glucose_rxn = "EX_glc__D_e"
    if target_rxn not in model.reactions:
        raise ValueError(f"Exchange reaction {target_rxn} not found in the model.")

    # ----------------------
    # Wild-type baseline
    # ----------------------
    model.objective = biomass_rxn
    wt_sol    = model.optimize()
    wt_flux   = wt_sol.fluxes
    wt_growth = wt_sol.objective_value

    # ----------------------
    # Theoretical max succinate
    # ----------------------
    succ_model = model.copy()
    succ_model.objective = target_rxn
    theoretical_succ = succ_model.optimize().objective_value

    # ----------------------
    # OptKnock module
    # ----------------------
    viability_threshold = 0.01 * wt_growth
    optknock_module = sd.SDModule(
        model=model,
        module_type=sd.OPTKNOCK,
        inner_objective=biomass_rxn,
        outer_objective=target_rxn,
        constraints=[f"{biomass_rxn} >= {viability_threshold}"]
    )

    # ----------------------
    # Run OptKnock
    # ----------------------
    max_solutions = 1
    solutions = sd.compute_strain_designs(
        model,
        sd_modules=[optknock_module],
        max_solutions=max_solutions,
        max_cost=4,
        solution_approach=sd.BEST,
        solver="gurobi"
    )

    # ----------------------
    # Loop over solutions
    # ----------------------
    ko_list     = []
    found_valid = False
    for i in range(max_solutions):
        try:
            reaction_sd = solutions.get_reaction_sd(i)
        except Exception:
            break
        candidate_ko = []
        skip = False
        for rd in reaction_sd:
            for rxn_id, val in rd.items():
                if is_knockout_value(val):
                    if rxn_id in excluded_rxns:
                        print(f"Solution {i}: Skipped because knockout {rxn_id} is in exclusion list.")
                        skip = True
                    candidate_ko.append(rxn_id)
        if skip:
            continue
        ko_list     = candidate_ko
        found_valid = True
        print(f"Selected solution {i} with knockouts: {ko_list}")
        break

    if not found_valid or not ko_list:
        print("No valid OptKnock solution found (all contain excluded reactions).")
        print(f"\n  Total time: {time.time() - t_start:.2f}s")
        return

    # ----------------------
    # Build mutant
    # ----------------------
    mutant = model.copy()
    for rxn_id in ko_list:
        if rxn_id in mutant.reactions:
            mutant.reactions.get_by_id(rxn_id).knock_out()

    # ----------------------
    # Mutant simulation
    # ----------------------
    mutant.objective = biomass_rxn
    mut_sol   = mutant.optimize()
    mut_flux  = mut_sol.fluxes
    mut_succ  = mut_flux.get(target_rxn, float("nan"))
    mut_glc   = abs(mut_flux.get(glucose_rxn, float("nan")))
    mut_yield = round(mut_succ / mut_glc, 6) if mut_glc > 0 else float("nan")

    # ----------------------
    # Prepare Excel output
    # ----------------------
    all_kos    = turn_off_reactions + ko_list
    wt_df      = wt_flux.reset_index()
    wt_df.columns = ["reaction_id", "wt_flux"]
    mut_df     = mut_flux.reset_index()
    mut_df.columns = ["reaction_id", "mutant_flux"]
    compare_df = wt_df.merge(mut_df, on="reaction_id", how="outer")
    compare_df["delta(mut-wt)"] = compare_df["mutant_flux"] - compare_df["wt_flux"]

    ko_df = pd.DataFrame({"knockout_reaction_id": all_kos})

    summary_df = pd.DataFrame({
        "metric": [
            "WT_growth",
            "WT_succinate",
            "Mut_growth",
            "Mut_succinate",
            "Mut_succ_yield",
            "theoretical_max_succinate",
            "num_KOs_total",
            "turn_off_reactions",
            "optknock_KOs",
            "total_time_s",
        ],
        "value": [
            wt_growth,
            wt_flux.get(target_rxn, float("nan")),
            mut_sol.objective_value,
            mut_succ,
            mut_yield,
            theoretical_succ,
            len(all_kos),
            ", ".join(turn_off_reactions),
            ", ".join(ko_list),
            round(time.time() - t_start, 2),
        ]
    })

    out_file = "OptKnock_iAF1260_succinate.xlsx"
    with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
        summary_df.to_excel(writer, sheet_name="summary",         index=False)
        ko_df.to_excel(writer,      sheet_name="knockouts",       index=False)
        compare_df.to_excel(writer, sheet_name="flux_comparison", index=False)

    print(f"OptKnock computation complete. Results written to: {out_file}")
    print(f"\n  Total time: {time.time() - t_start:.2f}s")

if __name__ == "__main__":
    main()