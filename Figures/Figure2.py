import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np

# -----------------------------
# Data
# -----------------------------
cases = ['Case 1', 'Case 2', 'Case 3']
x = np.arange(len(cases))
w = 0.22

time_rk_gurobi = [134.53, 328.56, 62364.87]
time_rk_glpk   = [157.73, 282.40, 69363.09]
time_ok_gurobi = [418.56, 812.11, 90245.27]

yield_rk_gurobi = [0.921, 1.325, 1.311]
yield_rk_glpk   = [0.039, 1.321, 0.952]
yield_ok_gurobi = [0.921, 1.325, 0.959]

# -----------------------------
# Colors
# -----------------------------
c_rk_g = '#1b6ca8'
c_rk_l = '#7fb3d3'
c_ok_g = '#e8722a'

# -----------------------------
# Axis style
# -----------------------------
def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.4)
    ax.spines['bottom'].set_linewidth(1.4)
    ax.tick_params(axis='both', width=1.4, labelsize=13)
    for lbl in ax.get_xticklabels() + ax.get_yticklabels():
        lbl.set_fontweight('bold')

# -----------------------------
# Figure layout
# -----------------------------
fig = plt.figure(figsize=(17, 7.5))
gs = gridspec.GridSpec(
    2, 2,
    width_ratios=[1, 1],
    height_ratios=[1, 2],
    wspace=0.35,
    hspace=0.12
)

# =====================================================
# Panel A (Broken Y-axis)
# =====================================================
ax_top = fig.add_subplot(gs[0, 0])
ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)

for ax in [ax_top, ax_bot]:
    ax.bar(x - w, time_rk_gurobi, w,
           color=c_rk_g, edgecolor='black', linewidth=0.8)
    ax.bar(x, time_rk_glpk, w,
           color=c_rk_l, edgecolor='black', linewidth=0.8)
    ax.bar(x + w, time_ok_gurobi, w,
           color=c_ok_g, edgecolor='black', linewidth=0.8)

# Broken axis limits
ax_bot.set_ylim(0, 1200)
ax_top.set_ylim(55000, 95000)

# Hide touching spines
ax_top.spines['bottom'].set_visible(False)
ax_bot.spines['top'].set_visible(False)

# Ticks
ax_top.tick_params(labelbottom=False, bottom=False)
ax_bot.set_xticks(x)
ax_bot.set_xticklabels(cases, fontsize=14, fontweight='bold')
ax_bot.set_ylabel('Time (s)', fontsize=15, fontweight='bold')
ax_top.set_title('A', loc='left', fontsize=17, fontweight='bold')

style_ax(ax_top)
style_ax(ax_bot)

# =====================================================
# Diagonal break marks (aspect-independent)
# =====================================================
d = 0.5  # slope proportion of the marker
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=10,
              linestyle="none", color='k', mec='k',
              mew=1.3, clip_on=False)
ax_top.plot([0, 1], [0, 0], transform=ax_top.transAxes, **kwargs)
ax_bot.plot([0, 1], [1, 1], transform=ax_bot.transAxes, **kwargs)

# =====================================================
# Panel B
# =====================================================
ax = fig.add_subplot(gs[:, 1])
ax.bar(x - w, yield_rk_gurobi, w,
       color=c_rk_g, edgecolor='black', linewidth=0.8)
ax.bar(x, yield_rk_glpk, w,
       color=c_rk_l, edgecolor='black', linewidth=0.8)
ax.bar(x + w, yield_ok_gurobi, w,
       color=c_ok_g, edgecolor='black', linewidth=0.8)

ax.set_ylabel('Yield (mol succinate / mol glucose)',
              fontsize=15, fontweight='bold')
ax.set_title('B', loc='left', fontsize=17, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(cases, fontsize=14, fontweight='bold')
ax.set_ylim(0, 1.4)
style_ax(ax)

# =====================================================
# Legend
# =====================================================
patches = [
    mpatches.Patch(facecolor=c_rk_g, edgecolor='black',
                   label='RapidKnock (Gurobi)'),
    mpatches.Patch(facecolor=c_rk_l, edgecolor='black',
                   label='RapidKnock (GLPK)'),
    mpatches.Patch(facecolor=c_ok_g, edgecolor='black',
                   label='OptKnock (Gurobi)')
]
fig.legend(handles=patches,
           loc='lower center',
           ncol=3,
           frameon=False,
           prop={'weight': 'bold', 'size': 13})

plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig('Figure2.png', dpi=300, bbox_inches='tight')
print("Saved Figure2.png")