import cobra
from cobra.flux_analysis import production_envelope
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

MODEL_FILE  = "iAF1260.xml"
BIOMASS_RXN = "BIOMASS_Ec_iAF1260_core_59p81M"
TARGET_RXN  = "EX_succ_e"
OXYGEN_RXN  = "EX_o2_e"
POINTS      = 30

ko_sets = {
    ("Case1", "OptKnock",   "Gurobi"): ["LDH_D", "ALCD2x"],
    ("Case1", "RapidKnock", "Gurobi"): ["LDH_D", "ALCD2x"],
    ("Case1", "RapidKnock", "GLPK"):   ["ATPS4rpp", "ALCD2x"],
    ("Case2", "OptKnock",   "Gurobi"): ["ACALD", "ATPS4rpp", "LDH_D", "ALCD2x"],
    ("Case2", "RapidKnock", "Gurobi"): ["PFL", "THD2pp", "LDH_D", "ACALD"],
    ("Case2", "RapidKnock", "GLPK"):   ["GLCptspp", "PYK", "TKT2", "ALCD2x"],
}

base_model = cobra.io.read_sbml_model(MODEL_FILE)
base_model.reactions.get_by_id(OXYGEN_RXN).lower_bound = 0

envelopes = {}
for key, kos in ko_sets.items():
    with base_model as model:
        for rxn_id in kos:
            model.reactions.get_by_id(rxn_id).knock_out()
        env = production_envelope(
            model, reactions=[BIOMASS_RXN], objective=TARGET_RXN, points=POINTS
        )
        envelopes[key] = env
        print(f"{key}: {len(env)} points")

def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.4)
    ax.spines['bottom'].set_linewidth(1.4)
    ax.tick_params(axis='both', width=1.4, labelsize=12)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight('bold')

c_ok = '#e8722a'
c_rk = '#1b6ca8'

panels = [
    ("A", "Case1", "Gurobi"),
    ("B", "Case1", "GLPK"),
    ("C", "Case2", "Gurobi"),
    ("D", "Case2", "GLPK"),
]

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

for ax, (label, case, rk_solver) in zip(axes, panels):
    ok_env = envelopes[(case, "OptKnock", "Gurobi")]
    rk_env = envelopes[(case, "RapidKnock", rk_solver)]

    # OptKnock drawn first, dashed + thicker, so it stays visible even when
    # RapidKnock's envelope fully overlaps or covers it
    ox = ok_env[BIOMASS_RXN]
    ax.plot(ox, ok_env["flux_maximum"], color=c_ok, linewidth=3.2,
            linestyle='--', label="OptKnock (Gurobi)", zorder=2)
    ax.plot(ox, ok_env["flux_minimum"], color=c_ok, linewidth=3.2,
            linestyle='--', zorder=2)
    ax.fill_between(ox, ok_env["flux_minimum"], ok_env["flux_maximum"],
                     color=c_ok, alpha=0.12, zorder=1)
    # close the OptKnock envelope at its max-biomass edge
    ax.plot([ox.iloc[-1], ox.iloc[-1]],
            [ok_env["flux_minimum"].iloc[-1], ok_env["flux_maximum"].iloc[-1]],
            color=c_ok, linewidth=3.2, linestyle='--', zorder=2)

    # RapidKnock drawn on top, solid + thinner, so both stay distinguishable
    rx = rk_env[BIOMASS_RXN]
    ax.plot(rx, rk_env["flux_maximum"], color=c_rk, linewidth=2.0,
            linestyle='-', label=f"RapidKnock ({rk_solver})", zorder=3)
    ax.plot(rx, rk_env["flux_minimum"], color=c_rk, linewidth=2.0,
            linestyle='-', zorder=3)
    ax.fill_between(rx, rk_env["flux_minimum"], rk_env["flux_maximum"],
                     color=c_rk, alpha=0.12, zorder=1)

    ax.set_xlabel('Biomass (1/h)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Succinate flux (mmol/gDW/h)', fontsize=13, fontweight='bold')
    ax.set_title(label, loc='left', fontsize=17, fontweight='bold')
    ax.legend(frameon=False, fontsize=10, loc='upper right')
    ax.margins(x=0.03, y=0.08)
    style_ax(ax)

plt.tight_layout()
plt.savefig('envelope_curves.png', dpi=300, bbox_inches='tight')
print("Saved envelope_curves.png")
