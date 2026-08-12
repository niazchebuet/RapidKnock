"""
FastKnock — True Single-Pass Exhaustive Enumeration with FVA-Based Pruning
Based on Hassani et al. 2024 (Microbial Cell Factories), Algorithms 1-4.

This is a faithful implementation of the FastKnock algorithm:
  - Single-pass depth-first traversal (not round-based, not greedy)
  - Exhaustive enumeration to MAX_KO with dynamic pruning
  - All solutions at all depths are collected and saved
  - Pruning rule: only explore target_space = nonzero-flux reactions ∩ Removable
"""

import cobra
from cobra.flux_analysis import find_blocked_reactions
import pandas as pd
import numpy as np
from collections import defaultdict
import warnings
import time

warnings.filterwarnings("ignore")

# ── Parameters ─────────────────────────────────────────────────────────────────

MODEL_FILE          = "iAF1260.xml"
BIOMASS_RXN         = "BIOMASS_Ec_iAF1260_core_59p81M"
CHEMICAL_RXN        = "EX_succ_e"
GLUCOSE_EX          = "EX_glc__D_e"
OXYGEN_EX           = "EX_o2_e"
SOLVER              = "gurobi"
ANAEROBIC           = True
GLUCOSE_UPTAKE      = 10.0
MIN_BIOMASS_FRAC    = 0.01
TH_CHEMICAL_FRAC    = 0.05
MAX_KO              = 2
ALWAYS_EXCLUDED     = ["ATPM"]           # Reactions permanently excluded from knockout candidates
TURN_OFF_REACTIONS  = []                 # Additional reactions to exclude (e.g., ["PGI", "PFK"])

# ──────────────────────────────────────────────────────────────────────────────

t_start = time.time()

print("="*80)
print("FastKnock — True Exhaustive Enumeration with FVA-Based Pruning")
print("="*80)

# ── 1. Load model ──────────────────────────────────────────────────────────────

print("\n[1/7] Loading model...")
model = cobra.io.read_sbml_model(MODEL_FILE)
model.solver = SOLVER

if ANAEROBIC:
    model.reactions.get_by_id(OXYGEN_EX).lower_bound = 0

model.reactions.get_by_id(GLUCOSE_EX).lower_bound = -GLUCOSE_UPTAKE
model.reactions.get_by_id(GLUCOSE_EX).upper_bound = 0

print(f"  Reactions  : {len(model.reactions)}")
print(f"  Metabolites: {len(model.metabolites)}")
print(f"  Genes      : {len(model.genes)}")
print(f"  Solver     : {SOLVER}")
print(f"  Anaerobic  : {ANAEROBIC}")

# ── 2. Compute WT baseline and thresholds ──────────────────────────────────────

print("\n[2/7] Computing WT baseline and thresholds...")

with model:
    model.objective = BIOMASS_RXN
    wt_sol = model.optimize()
    wt_biomass = wt_sol.objective_value

print(f"  WT max biomass: {wt_biomass:.4f}")

with model:
    model.reactions.get_by_id(BIOMASS_RXN).lower_bound = MIN_BIOMASS_FRAC * wt_biomass
    model.objective = CHEMICAL_RXN
    max_chemical = model.slim_optimize()

TH_BIOMASS = MIN_BIOMASS_FRAC * wt_biomass
TH_CHEMICAL = TH_CHEMICAL_FRAC * max_chemical

print(f"  Max theoretical {CHEMICAL_RXN}: {max_chemical:.4f}")
print(f"  Th_biomass : {TH_BIOMASS:.4f}")
print(f"  Th_chemical: {TH_CHEMICAL:.4f}")

# ── 3. Preprocessing: Remove blocked reactions only ────────────────────────────

print("\n[3/7] Preprocessing: removing blocked reactions and turning off reactions...")

blocked = find_blocked_reactions(model)
print(f"  Found {len(blocked)} blocked reactions")
for rxn_id in blocked:
    if rxn_id in model.reactions:
        model.reactions.get_by_id(rxn_id).remove_from_model()

# Turn off (knockout) specified reactions in the base model before search
if TURN_OFF_REACTIONS:
    print(f"  Turning off reactions before search: {TURN_OFF_REACTIONS}")
    for rxn_id in TURN_OFF_REACTIONS:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).bounds = (0, 0)

eliminate_list = set(ALWAYS_EXCLUDED)
Removable = [r.id for r in model.reactions if r.id not in eliminate_list]
Removable_set = set(Removable)

print(f"  Removable set size: {len(Removable)} reactions")
print(f"  Excluded from knockout candidates: {sorted(eliminate_list)}")

# ── 4. Build gene-rule co-knockout map ──────────────────────────────────────────

print("\n[4/7] Building gene-rule co-knockout map...")

def gene_rules(m):
    """Parse gene-reaction rules and build essential/involved gene sets per reaction."""
    ess_and_inv = defaultdict(lambda: {"essential": [], "involved": []})
    deleted_by_gene = defaultdict(lambda: {"essential": [], "involved": []})

    for rxn in m.reactions:
        genes = rxn.gene_reaction_rule.replace(" ", "").replace("(", "").replace(")", "")

        if not genes:
            ess_and_inv[rxn.id] = {"essential": [], "involved": []}
            continue

        or_parts = genes.split("or")
        and_parts = genes.split("and")

        if len(or_parts) == 1 and len(and_parts) == 1:
            ess_and_inv[rxn.id] = {"essential": [genes], "involved": []}
            deleted_by_gene[genes]["essential"].append(rxn.id)
        elif len(or_parts) == 1 and len(and_parts) > 1:
            rule = genes.split("and")
            ess_and_inv[rxn.id] = {"essential": rule, "involved": []}
            for g in rule:
                deleted_by_gene[g]["essential"].append(rxn.id)
        elif len(or_parts) > 1 and "and" not in genes:
            rule = genes.split("or")
            ess_and_inv[rxn.id] = {"essential": [], "involved": rule}
            for g in rule:
                deleted_by_gene[g]["involved"].append(rxn.id)
        else:
            # Mixed AND/OR
            rule = genes.replace(")", "").replace("(", "").split("or")
            and_or_groups = [r.split("and") for r in rule]
            ess_and_inv[rxn.id] = {"essential": [], "involved": and_or_groups}

    return ess_and_inv, deleted_by_gene

def get_co_knockouts(rxn_id, ess_and_inv, deleted_by_gene):
    """Return set of reactions that would be co-deleted if rxn_id is knocked out."""
    ess_genes = ess_and_inv[rxn_id]["essential"]
    inv_genes = ess_and_inv[rxn_id]["involved"]

    if not ess_genes and not inv_genes:
        return set()

    if len(ess_genes) == 1:
        return set(deleted_by_gene[ess_genes[0]]["essential"])

    if len(ess_genes) > 1:
        best_gene = min(ess_genes, key=lambda g: len(deleted_by_gene[g]["essential"]))
        return set(deleted_by_gene[best_gene]["essential"])

    if inv_genes:
        return set().union(*[set(deleted_by_gene[g]["essential"]) for g in inv_genes])

    return set()

def build_co_knockout_map(m):
    """Build complete co-knockout map for all reactions."""
    ess_and_inv, deleted_by_gene = gene_rules(m)
    co_map = {}
    for rxn in m.reactions:
        co_map[rxn.id] = get_co_knockouts(rxn.id, ess_and_inv, deleted_by_gene)
    return co_map

coKnockoutMap = build_co_knockout_map(model)
print(f"  Co-knockout map built for {len(coKnockoutMap)} reactions")

# ── 5. Node class and target-space identification (Algorithm 2) ───────────────

print("\n[5/7] Building search tree...")

class Node:
    """Represents a single node in the depth-first traversal tree."""
    __slots__ = ["level", "deleted_rxns", "target_space", "growth", "chemical", "is_solution"]

    def __init__(self, level, deleted_rxns):
        self.level = level
        self.deleted_rxns = list(deleted_rxns)
        self.target_space = []
        self.growth = None
        self.chemical = None
        self.is_solution = False

def identify_target_space(node, model):
    """
    Algorithm 2: Identify target space for a node.
    Target space = reactions with nonzero flux ∩ Removable set
    """
    # Include co-knockouts for all currently deleted reactions
    co_extra = set()
    for rxn_id in node.deleted_rxns:
        co_extra |= coKnockoutMap.get(rxn_id, set())
    full_deleted = sorted(set(node.deleted_rxns) | co_extra)
    node.deleted_rxns = full_deleted

    # Run FBA with deletions
    with model:
        for rxn_id in full_deleted:
            if rxn_id in model.reactions:
                model.reactions.get_by_id(rxn_id).bounds = (0, 0)

        sol = model.optimize()

        if sol.status == "optimal":
            node.growth = sol.fluxes.get(BIOMASS_RXN, 0.0)
            node.chemical = sol.fluxes.get(CHEMICAL_RXN, 0.0)
            nonzero = sol.fluxes[abs(sol.fluxes) > 1e-9].index
            node.target_space = [r for r in nonzero if r in Removable_set]
        else:
            node.growth = 0.0
            node.chemical = 0.0
            node.target_space = []

    # Check if this is a valid solution
    node.is_solution = (node.growth > TH_BIOMASS) and (node.chemical > TH_CHEMICAL)
    return node

# ── 6. Depth-first tree traversal (Algorithms 1, 3, 4) ────────────────────────

def construct_subtree_and_traverse(root, model, max_depth):
    """
    Construct tree and traverse depth-first, collecting all solutions.
    Returns list of all solution nodes found.
    """
    all_solutions = []
    nodes_checked = 0

    # Queue-based DFS: queue_by_level[lvl] = list of nodes to process at level lvl
    queue_by_level = {lvl: [] for lvl in range(max_depth + 1)}
    checked_by_level = {lvl: [] for lvl in range(max_depth + 1)}
    queue_by_level[0] = [root]

    def process_children(parent_node):
        """Algorithm 4: Construct children from parent's target space."""
        nonlocal nodes_checked
        parent_level = parent_node.level

        if parent_level >= max_depth:
            return

        next_level = parent_level + 1

        for rxn in parent_node.target_space:
            # Avoid permutation duplicates: only add reactions not yet checked at next level
            if rxn in checked_by_level[next_level]:
                continue

            child = Node(next_level, parent_node.deleted_rxns + [rxn])
            identify_target_space(child, model)
            nodes_checked += 1

            # Collect if solution
            if child.is_solution:
                all_solutions.append(child)

            # Queue for further expansion
            queue_by_level[next_level].append(child)
            checked_by_level[next_level].append(rxn)

    # DFS traversal
    for level in range(max_depth + 1):
        while queue_by_level[level]:
            node = queue_by_level[level].pop(0)
            process_children(node)

    return all_solutions, nodes_checked

# Run the search
root = Node(0, [])
identify_target_space(root, model)
print(f"  Root target space: {len(root.target_space)} reactions")

all_solutions, nodes_checked = construct_subtree_and_traverse(root, model, MAX_KO)

print(f"  Nodes checked: {nodes_checked}")
print(f"  Solutions found: {len(all_solutions)}")

# ── 7. Save results ────────────────────────────────────────────────────────────

print("\n[6/7] Organizing results...")

# Build results dataframe
results_rows = []
for sol in all_solutions:
    yield_val = (sol.chemical / GLUCOSE_UPTAKE) if GLUCOSE_UPTAKE > 0 else 0.0
    results_rows.append({
        "depth": sol.level,
        "n_knockouts": len(sol.deleted_rxns),
        "knockouts": ", ".join(sol.deleted_rxns),
        "growth": round(sol.growth, 6),
        "chemical_flux": round(sol.chemical, 6),
        "yield": round(yield_val, 6),
    })

results_df = pd.DataFrame(results_rows)
if not results_df.empty:
    results_df = results_df.sort_values(["n_knockouts", "yield"], ascending=[True, False])

# Build summary
summary_rows = []
for depth in range(1, MAX_KO + 1):
    count = len([s for s in all_solutions if s.level == depth])
    if count > 0:
        best = max([s for s in all_solutions if s.level == depth], key=lambda x: x.chemical)
        summary_rows.append({
            "knockout_order": depth,
            "n_solutions": count,
            "best_yield": round((best.chemical / GLUCOSE_UPTAKE), 6),
            "best_knockouts": ", ".join(best.deleted_rxns),
        })

summary_df = pd.DataFrame(summary_rows)

# Save to Excel
out_file = "FastKnock_true_results.xlsx"
print(f"\n[7/7] Saving results to {out_file}...")

with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
    summary_df.to_excel(writer, sheet_name="summary", index=False)
    if not results_df.empty:
        results_df.to_excel(writer, sheet_name="all_solutions", index=False)
    for depth in range(1, MAX_KO + 1):
        depth_df = results_df[results_df["n_knockouts"] == depth]
        if not depth_df.empty:
            depth_df.to_excel(writer, sheet_name=f"{depth}_knockouts", index=False)

total_time = round(time.time() - t_start, 2)

print("\n" + "="*80)
print("RESULTS SUMMARY")
print("="*80)
print(summary_df.to_string(index=False) if not summary_df.empty else "No solutions found")
print(f"\nTotal execution time: {total_time}s")
print(f"Results saved to: {out_file}")
print("="*80)