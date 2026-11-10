"""Main and supplementary figures from the analysis-module outputs."""
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, Rectangle, Patch
from matplotlib.lines import Line2D
OUT = Path('output/figures')
TABLES = Path('output/tables')
C = {'navy': '#24485C', 'teal': '#2C7C87', 'copper': '#C07852', 'purple': '#8C7AA7', 'sage': '#8B9D77', 'gray': '#5F6570', 'pale': '#E8EDF0', 'ink': '#1E2930'}

def setup():
    for path in Path('assets/fonts').glob('*'):
        if path.suffix.lower() == '.ttf':
            font_manager.fontManager.addfont(str(path))
    font = font_manager.findfont(font_manager.FontProperties(family='Times New Roman'), fallback_to_default=False)
    if font_manager.FontProperties(fname=font).get_name() != 'Times New Roman':
        raise ValueError('Times New Roman is required')
    plt.rcParams.update({'font.family': 'Times New Roman', 'font.size': 9.5, 'axes.titlesize': 12, 'axes.labelsize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.65, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 3, 'ytick.major.size': 3, 'axes.edgecolor': C['gray'], 'axes.labelcolor': C['ink'], 'text.color': C['ink'], 'xtick.color': C['ink'], 'ytick.color': C['ink'], 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white', 'figure.facecolor': 'white', 'axes.unicode_minus': True, 'mathtext.fontset': 'custom', 'mathtext.rm': 'Times New Roman', 'mathtext.it': 'Times New Roman:italic', 'mathtext.bf': 'Times New Roman:bold'})
    return font

def label(ax, letter, title):
    ax.set_title(title, loc='left', pad=16, fontweight='bold', fontsize=12)
    ax.text(-0.13, 1.07, letter, transform=ax.transAxes, fontsize=14, fontweight='bold', va='bottom')

def save(fig, stem):
    OUT.mkdir(parents=True, exist_ok=True)
    for extension in ('pdf', 'svg'):
        fig.savefig(OUT / f'{stem}.{extension}', bbox_inches='tight', pad_inches=0.18)
    fig.savefig(OUT / f'{stem}.tiff', dpi=600, bbox_inches='tight', pad_inches=0.18, pil_kwargs={'compression': 'tiff_lzw'})
    fig.savefig(OUT / f'{stem}_preview.png', dpi=150, bbox_inches='tight', pad_inches=0.18)
    plt.close(fig)

def read_airway(context, filename):
    return pd.read_csv(TABLES / 'airway' / context / filename, sep='\t')

def figure_1():
    D = json.loads(Path('data/figure_metadata/evidence_inventory.json').read_text(encoding='utf-8'))
    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(3, 2, left=0.075, right=0.975, top=0.89, bottom=0.11, hspace=0.66, wspace=0.35, height_ratios=[1, 1.02, 1])
    axs = [fig.add_subplot(gs[i, j]) for i in range(3) for j in range(2)]

    def title(ax, letter, text):
        ax.set_title(text, loc='left', pad=12, fontsize=12, fontweight='bold')
        ax.text(-0.105, 1.055, letter, fontsize=15, fontweight='bold', va='bottom', transform=ax.transAxes)

    def clean(ax):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
    ax = axs[0]
    clean(ax)
    title(ax, 'A', 'The hypothesis requires links between layers')
    nodes = [(0.02, 0.69, 'Clinical phenotype', 'Explicit cough status\nand duration', C['navy']), (0.58, 0.69, 'Molecular context', 'Blood proteins and\nnasal transcriptomes', C['teal']), (0.58, 0.16, 'Local physiology', 'Laryngeal motor activity;\nsensory pathway unmeasured', C['purple']), (0.02, 0.16, 'Patient benefit', 'Daily cough and quality of life;\nbenefit studied in separate cohorts', C['copper'])]
    for x, y, heading, detail, col in nodes:
        ax.plot([x, x], [y - 0.015, y + 0.23], color=col, lw=3)
        ax.text(x + 0.025, y + 0.22, heading, color=col, fontsize=11.4, fontweight='bold', va='top')
        ax.text(x + 0.025, y + 0.105, detail, fontsize=10.1, va='top', linespacing=1.35)
    for start, end in [((0.43, 0.8), (0.555, 0.8)), ((0.78, 0.65), (0.78, 0.45)), ((0.55, 0.27), (0.43, 0.27))]:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=12, linestyle=(0, (3, 3)), color=C['gray'], lw=1.1))
    ax.text(0.4, 0.5, 'No matched analysis of\nmolecular context, physiology\nand treatment response', ha='center', va='center', fontsize=9.7, fontstyle='italic', color=C['gray'])
    ax.plot([0.02, 0.12], [0.015, 0.015], ls=(0, (3, 3)), color=C['gray'], lw=1.1)
    ax.text(0.14, 0.015, 'Dashed arrows: conceptual links, not estimated causal effects', va='center', fontsize=9)
    ax = axs[1]
    title(ax, 'B', 'Participant numbers depend on the analysis')
    rows = D['counts']
    yy = np.arange(len(rows))[::-1]
    layer_colors = {'blood': C['teal'], 'airway': C['navy'], 'physiology': C['purple'], 'trial': C['copper']}
    ax.barh(yy, [r['n'] for r in rows], height=0.6, color=[layer_colors[r['layer']] for r in rows], alpha=0.95)
    for y, r in zip(yy, rows):
        ax.text(r['n'] + 2, y, str(r['n']), va='center', fontsize=10, fontweight='bold')
    ax.set_yticks(yy, [r['label'] for r in rows])
    ax.set_xlim(0, 118)
    ax.set_ylim(-0.6, 8.6)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel('Participants or donors (not libraries)', labelpad=3)
    ax.tick_params(axis='y', length=0, pad=5)
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x', color=C['pale'], linewidth=0.6)
    ax.set_axisbelow(True)
    ax.text(0, -0.27, 'INCOV + MY-LC: 203 eligible → 202 with protein data.\nLongitudinal/radiology analyses share a cohort; no pooled total N.', transform=ax.transAxes, fontsize=9.2, va='top', linespacing=1.3)
    ax = axs[2]
    title(ax, 'C', 'Sampling time is not cough duration')
    clean(ax)
    ax.text(0, 0.98, 'Observed postinfection assessment ranges', fontsize=10.5, fontweight='bold', va='top')
    chart = ax.inset_axes([0.12, 0.52, 0.85, 0.35])
    chart.set_xlim(0, 620)
    chart.set_ylim(-0.55, 1.55)
    for y, r, col in zip([1, 0], D['observed_timing'], [C['teal'], C['navy']]):
        chart.plot([r['min_days'], r['max_days']], [y, y], lw=5, color=col, solid_capstyle='butt')
        chart.scatter([r['min_days'], r['max_days']], [y, y], s=23, color=col, zorder=3)
        label = f"{r['min_days']:.2f}–{r['max_days']:.2f} d" if r['cohort'] == 'INCOV' else '104–578 d'
        chart.text((r['min_days'] + r['max_days']) / 2, y + 0.2, label, ha='center', fontsize=9.5)
    chart.set_yticks([1, 0], ['INCOV', 'MY-LC'])
    chart.tick_params(axis='y', length=0, pad=5)
    chart.set_xticks([0, 200, 400, 600])
    chart.set_xlabel('Days after recorded symptom onset', fontsize=9, labelpad=2)
    chart.spines['left'].set_visible(False)
    chart.spines['bottom'].set_color(C['gray'])
    ax.plot([0, 1], [0.34, 0.34], color=C['pale'], lw=1)
    ax.text(0, 0.28, 'Study entry criteria (different constructs)', fontsize=10.5, fontweight='bold', va='top')
    ax.text(0, 0.145, 'LEMG', fontsize=10, color=C['purple'], fontweight='bold', va='center')
    ax.text(0.21, 0.145, 'Cough present ≥12 weeks after COVID-19', fontsize=9.7, va='center')
    ax.text(0, 0.025, 'Selected trial', fontsize=10, color=C['copper'], fontweight='bold', va='center')
    ax.text(0.21, 0.025, 'Daily dry cough for ≥6 months', fontsize=9.7, va='center')
    ax.text(0, -0.105, 'Neither sampling range establishes continuous cough in each participant.', fontsize=9.2, va='top', color=C['gray'])
    ax = axs[3]
    title(ax, 'D', 'Tissue and endpoint coverage are uneven')
    arr = np.array(D['availability'])
    nr, nc = arr.shape
    for i in range(nr):
        for j in range(nc):
            ax.add_patch(Rectangle((j - 0.47, i - 0.43), 0.94, 0.86, facecolor='#F5F7F8', edgecolor='white', lw=0.5))
            if arr[i, j] == 1:
                ax.scatter(j, i, s=52, color=C['teal'], zorder=3)
            elif arr[i, j] == 2:
                ax.scatter(j, i, s=62, facecolors='none', edgecolors=C['copper'], linewidths=1.5, zorder=3)
            else:
                ax.text(j, i, '–', ha='center', va='center', fontsize=11, color='#A9AFB3')
    ax.set_xlim(-0.6, 5.6)
    ax.set_ylim(8.55, -2.0)
    ax.set_yticks(np.arange(nr), D['availability_rows'])
    ax.set_xticks([])
    for j, label in enumerate(D['availability_columns']):
        ax.text(j, -1.3, label, ha='center', va='center', fontsize=8.8)
    ax.tick_params(axis='both', length=0, pad=5, labelsize=8.8)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(0, -0.12, '● available   ○ incomplete cough coding   – not in analysed inputs', fontsize=9.2, transform=ax.transAxes, va='top')
    ax.text(0, -0.23, 'Blood RNA: “not listed” does not establish absence of cough.', fontsize=9.2, transform=ax.transAxes, va='top', color=C['gray'])
    ax = axs[4]
    clean(ax)
    title(ax, 'E', 'Accessions and analyses are not independent units')

    def source(x, y, text, w=0.36):
        ax.add_patch(Rectangle((x - w / 2, y - 0.085), w, 0.17, facecolor='#F1F4F5', edgecolor=C['pale'], lw=0.8))
        ax.text(x, y, text, ha='center', va='center', fontsize=9.5)

    def analysis(y, text, col):
        ax.scatter(0.57, y, s=38, color=col, zorder=4)
        ax.text(0.6, y, text, ha='left', va='center', color=col, fontsize=10.5, fontweight='bold', zorder=5)
    source(0.2, 0.84, 'GSE286429')
    analysis(0.84, 'Culture contrasts', C['navy'])
    ax.plot([0.38, 0.55], [0.84, 0.84], color=C['navy'], lw=1.1)
    ax.text(0.76, 0.7, 'same study / donors', ha='center', fontsize=9, color=C['gray'])
    source(0.2, 0.49, 'GSE261777')
    ax.plot([0.38, 0.49], [0.49, 0.49], color=C['navy'], lw=1.1)
    ax.plot([0.49, 0.49], [0.41, 0.57], color=C['navy'], lw=1.1)
    for y, text in [(0.57, 'Longitudinal'), (0.41, 'Radiology')]:
        analysis(y, text, C['navy'])
        ax.plot([0.49, 0.55], [y, y], color=C['navy'], lw=1.1)
    ax.text(0.76, 0.29, 'same source cohort', ha='center', fontsize=9, color=C['gray'])
    source(0.2, 0.14, 'GSE232623\nGSE232627')
    analysis(0.14, 'Pseudobulk', C['teal'])
    ax.plot([0.38, 0.55], [0.14, 0.14], color=C['teal'], lw=1.1)
    ax.plot([0.2, 0.2], [0.58, 0.75], ls=(0, (3, 3)), lw=1.1, color=C['copper'])
    ax.text(0.02, 0.665, 'overlap\nunresolved', ha='left', va='center', fontsize=8.9, color=C['copper'])
    ax.text(0, -0.08, 'Solid: shared source data. Dashed: culture–clinical overlap unresolved.', fontsize=9.2, va='top', color=C['gray'])
    ax = axs[5]
    title(ax, 'F', 'Multiplicity is controlled within defined families')
    labels = ['Context analyses', 'Blood primary', 'INCOV exploratory', 'MY-LC exploratory']
    values = [50, 12, 448, 138]
    ys = [3, 2, 1, 0]
    colors = [C['navy'], C['teal'], C['purple'], C['sage']]
    ax.barh(ys, values, height=0.52, color=colors)
    for boundary in [10, 20, 30, 40]:
        ax.plot([boundary, boundary], [2.74, 3.26], color='white', lw=1)
    for y, n, family in zip(ys, values, ['5 × BH(10)', 'BH(12)', 'BH(448)', 'BH(138)']):
        ax.text(n + 8, y, f'{n}  |  {family}', va='center', fontsize=10)
    ax.set_xlim(0, 610)
    ax.set_ylim(-0.6, 3.6)
    ax.set_yticks(ys, labels)
    ax.set_xticks([0, 100, 200, 300, 400, 500])
    ax.set_xlabel('Number of planned or screened tests', labelpad=5)
    ax.tick_params(axis='y', length=0, pad=5)
    ax.spines['left'].set_visible(False)
    ax.grid(axis='x', color=C['pale'], lw=0.6)
    ax.set_axisbelow(True)
    ax.text(0, -0.27, 'Context: five analyses × ten programmes, with separate BH correction.\nExploratory: 586 tests across two families; no pooled 586-test q value.', transform=ax.transAxes, fontsize=9.2, va='top', linespacing=1.3)
    save(fig, 'Figure1')

def figure_2():
    inp = TABLES / 'proteomics'
    elig = pd.read_csv('data/processed/proteomics/eligibility.csv')
    primary = pd.read_csv(inp / 'primary_panel_results.csv')
    exploratory = pd.read_csv(inp / 'exploratory_assay_results.csv')
    prec = pd.read_csv(inp / 'primary_precision.csv')
    assert len(primary) == 12 and len(exploratory) == 586
    eligible = elig[elig.eligible.astype(str).str.lower().eq('true')].copy()
    assert eligible.groupby('cohort').size().to_dict() == {'INCOV': 105, 'MYLC': 98}
    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(2, 3, left=0.072, right=0.975, bottom=0.135, top=0.855, wspace=0.38, hspace=0.66)
    axs = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]
    a = axs[0]
    label(a, 'A', 'From public records to analysis')
    a.axis('off')
    a.set_xlim(-0.45, 3.45)
    a.set_ylim(-0.65, 2.55)
    for x, t in enumerate(['Public\nrecords', 'Clinical\nsubset', 'Eligible\nparticipants', 'Primary\nmodels']):
        a.text(x, 2.3, t, ha='center', va='top', fontsize=9.5, color=C['gray'])
    for y, cohort, nums, color, sub in [(1.4, 'INCOV', [127, 106, 105, 105], C['navy'], 'Cough known'), (0.15, 'MYLC', [185, 101, 98, 97], C['teal'], 'Long COVID')]:
        a.text(-0.38, y + 0.35, cohort, color=color, fontweight='bold', fontsize=11)
        for x, n in enumerate(nums):
            a.text(x, y, str(n), ha='center', va='center', fontsize=20, fontweight='bold', color=color)
            if x < 3:
                a.annotate('', xy=(x + 0.78, y), xytext=(x + 0.22, y), arrowprops={'arrowstyle': '->', 'lw': 0.8, 'color': C['gray']})
        a.text(1, y - 0.31, sub, ha='center', fontsize=9)
        a.text(3, y - 0.31, '26 / 79' if cohort == 'INCOV' else '24 / 73', ha='center', fontsize=9)
    a.text(-0.38, -0.48, 'Model counts below: cough / no cough.\nMYLC: one eligible participant lacked all protein outcomes.', fontsize=9, color=C['gray'], va='top')
    a = axs[1]
    label(a, 'B', 'Timing differs substantially by cohort')
    rng = np.random.default_rng(20482)
    for i, (cohort, cough) in enumerate([('INCOV', 0), ('INCOV', 1), ('MYLC', 0), ('MYLC', 1)]):
        d = eligible[(eligible.cohort == cohort) & (eligible.cough == cough)].days.to_numpy()
        color = C['navy'] if cohort == 'INCOV' else C['teal']
        a.scatter(i + rng.uniform(-0.2, 0.2, len(d)), d, s=13, color=color, alpha=0.57, linewidth=0, zorder=2)
        q1, med, q3 = np.quantile(d, [0.25, 0.5, 0.75])
        a.vlines(i, q1, q3, color=color, lw=3, zorder=3)
        a.hlines(med, i - 0.23, i + 0.23, color='white', lw=3.4, zorder=4)
        a.hlines(med, i - 0.23, i + 0.23, color=color, lw=1.5, zorder=5)
        a.text(i, -0.13, f'n = {len(d)}', ha='center', transform=a.get_xaxis_transform(), fontsize=9)
    a.set_xticks(range(4), ['No cough', 'Cough', 'No cough', 'Cough'])
    a.set_ylim(0, 650)
    a.set_ylabel('Days since reported illness anchor')
    a.set_yticks([0, 200, 400, 600])
    a.set_xlim(-0.5, 3.5)
    a.text(0.125, -0.24, 'INCOV', transform=a.transAxes, color=C['navy'], ha='center', fontweight='bold')
    a.text(0.79, -0.24, 'MYLC', transform=a.transAxes, color=C['teal'], ha='center', fontweight='bold')
    a.grid(axis='y', color=C['pale'], lw=0.65, zorder=0)
    a = axs[2]
    label(a, 'C', 'No primary contrast survives correction')
    proteins = ['IFNG', 'TNF', 'IL10', 'CXCL9', 'CXCL10', 'CXCL11']
    for ci, cohort in enumerate(['INCOV', 'MYLC']):
        for j, p in enumerate(proteins):
            row = primary[(primary.cohort == cohort) & (primary.protein == p)].iloc[0]
            y = 5 - j + (0.13 if ci == 0 else -0.13)
            color = C['navy'] if ci == 0 else C['teal']
            a.errorbar(row.beta, y, xerr=[[row.beta - row.ci_low], [row.ci_high - row.beta]], fmt='o' if ci == 0 else 's', markersize=4.5, lw=1.15, capsize=2, color=color)
    a.axvline(0, color=C['gray'], lw=0.8, ls=(0, (3, 3)))
    a.set_yticks(range(6), proteins[::-1])
    a.set_ylim(-0.6, 5.85)
    a.set_xlim(-0.65, 1.15)
    a.set_xlabel('Adjusted difference / outcome SD')
    a.legend(handles=[Line2D([0], [0], marker='o', color=C['navy'], label='INCOV', lw=1), Line2D([0], [0], marker='s', color=C['teal'], label='MYLC', lw=1)], loc='upper center', bbox_to_anchor=(0.5, 1.08), ncol=2, frameon=False, handlelength=1.5, columnspacing=1)
    a.text(0.5, -0.23, 'HC3 95% CIs; minimum BH q = 0.300', transform=a.transAxes, ha='center', fontsize=9, color=C['gray'])
    for a, cohort, letter, n, color in [(axs[3], 'INCOV', 'D', 448, C['navy']), (axs[4], 'MYLC', 'E', 138, C['teal'])]:
        d = exploratory[exploratory.cohort == cohort]
        assert len(d) == n
        label(a, letter, f'{cohort}: all {n} exploratory assays')
        a.scatter(d.beta, -np.log10(d.p), s=16, c=color, alpha=0.52, linewidths=0.25, edgecolors='white')
        a.axhline(-np.log10(0.05), lw=0.8, color=C['gray'], ls=(0, (3, 3)))
        a.axvline(0, lw=0.6, color=C['gray'])
        a.text(0.97, -np.log10(0.05) + 0.08, 'Nominal P = 0.05', transform=a.get_yaxis_transform(), ha='right', fontsize=9, color=C['gray'])
        a.set_xlabel('Adjusted difference / outcome SD')
        a.set_ylabel('−log10 (nominal P)')
        a.set_ylim(0, max(-np.log10(d.p).max() + 0.35, 2))
        a.set_xlim(-max(abs(d.beta)) - 0.1, max(abs(d.beta)) + 0.1)
        a.text(0.03, 0.96, f'BH q < 0.05: 0 / {n}\nMinimum q = {d.q_exploratory.min():.3f}', transform=a.transAxes, va='top', fontsize=10, bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.88, 'pad': 3})
    a = axs[5]
    label(a, 'F', 'Non-significance leaves broad uncertainty')
    bounds = [0.2, 0.5, 0.8]
    nom = []
    fam = []
    for b in bounds:
        nom.append(int(((prec.ci_low <= b) & (prec.ci_high >= b)).sum()))
        fam.append(int(((prec.family95_low <= b) & (prec.family95_high >= b)).sum()))
    x = np.arange(3)
    a.bar(x - 0.18, nom, width=0.32, color=C['sage'], label='Nominal 95% CI')
    a.bar(x + 0.18, fam, width=0.32, color=C['purple'], label='Family-wise 95% CI')
    for xi, v in zip(x - 0.18, nom):
        a.text(xi, v + 0.23, str(v), ha='center', fontsize=10)
    for xi, v in zip(x + 0.18, fam):
        a.text(xi, v + 0.23, str(v), ha='center', fontsize=10)
    a.set_xticks(x, ['+0.2 SD', '+0.5 SD', '+0.8 SD'])
    a.set_ylim(0, 15.5)
    a.set_yticks([0, 4, 8, 12])
    a.set_ylabel('Primary CIs containing the reference value')
    a.set_xlabel('Illustrative positive effect reference')
    a.legend(loc='upper left', frameon=False, fontsize=9)
    a.text(0.5, -0.2, 'Reference values are not clinical thresholds.', transform=a.transAxes, ha='center', fontsize=9, color=C['gray'])
    save(fig, 'Figure2')

def figure_3():
    COL = ['#24485C', '#2C7C87', '#C07852', '#8C7AA7', '#8B9D77', '#5F6570']
    GRID = '#D9DFDF'
    INK = '#24353D'
    CMAP = LinearSegmentedColormap.from_list('airway_q', ['#FAFBFB', '#D8E6E6', COL[1], COL[0]])
    rows = []
    for code, context in [('B', 'culture'), ('D', 'longitudinal'), ('E', 'radiology'), ('F', 'pseudobulk'), ('G', 'blood')]:
        table = read_airway(context, 'primary_ROAST.tsv')
        table['analysis'] = code
        rows.append(table)
    primary = pd.concat(rows, ignore_index=True)
    changes = read_airway('longitudinal', 'program_donor_changes.tsv')
    summary = read_airway('longitudinal', 'program_donor_summary.tsv').set_index('gene_set')
    loo = read_airway('longitudinal', 'all_LOO_ROAST.tsv')
    es = pd.concat([read_airway('radiology', f'{name}_ROAST.tsv') for name in ['primary', 'composition37', 'RMSE34_composition']], ignore_index=True)
    es = es.loc[es.gene_set.eq('HALLMARK_TGF_BETA_SIGNALING')].reset_index(drop=True)
    donors = read_airway('pseudobulk', 'donor_manifest.tsv')
    programs = primary.loc[primary.analysis == 'B', 'gene_set'].tolist()
    labels = ['Cilium movement', 'Motile-cilium assembly', 'Tight-junction assembly', 'Epithelial differentiation', 'Mucus secretion', 'TNF-α / NF-κB', 'Interferon-α response', 'Interferon-γ response', 'TGF-β signalling', 'EMT']
    label = dict(zip(programs, labels))
    assert len(programs) == 10 and len(primary) == 50
    assert primary.groupby('analysis').size().eq(10).all()
    counts = {k: int((v.FDR < 0.05).sum()) for k, v in primary.groupby('analysis')}
    assert counts == {'B': 0, 'D': 4, 'E': 0, 'F': 0, 'G': 0}
    assert changes.groupby('gene_set').size().eq(6).all()
    assert loo.groupby('gene_set').size().eq(6).all()
    assert (summary.signflip_FDR > 0.05).all()
    assert donors.groupby('Status').size().to_dict() == {'Healthy': 6, 'RPRA': 5}
    assert np.allclose(es['FDR.Mixed'], [0.041, 0.1725, 0.3145])
    qmat = primary.pivot(index='gene_set', columns='analysis', values='FDR').loc[programs, ['B', 'D', 'E', 'F', 'G']]
    positive = primary[(primary.analysis == 'D') & (primary.FDR < 0.05)].gene_set.tolist()
    assert len(positive) == 4

    def head(ax, letter, title, subtitle=None):
        ax.set_title(title, loc='left', pad=20 if subtitle else 10, fontweight='bold', color=INK)
        ax.text(-0.16, 1.13 if subtitle else 1.075, letter, transform=ax.transAxes, fontsize=14, fontweight='bold', ha='left', va='bottom', color=INK)
        if subtitle:
            ax.text(0, 1.035, subtitle, transform=ax.transAxes, fontsize=9.2, color=COL[5], va='bottom')

    def polish(ax, xgrid=True):
        ax.tick_params(axis='both', length=3, pad=3, color=COL[5])
        ax.spines['left'].set_color(COL[5])
        ax.spines['bottom'].set_color(COL[5])
        if xgrid:
            ax.grid(axis='x', color=GRID, linewidth=0.45, zorder=0)

    def heat(ax, q, rows, cols, fontsize=9, vmax=2.35):
        z = -np.log10(q)
        im = ax.imshow(z, cmap=CMAP, vmin=0, vmax=vmax, aspect='auto', interpolation='nearest')
        ax.set_xticks(range(len(cols)), cols)
        ax.set_yticks(range(len(rows)), rows)
        ax.tick_params(length=0, pad=5)
        for s in ax.spines.values():
            s.set_visible(False)
        for i in range(q.shape[0]):
            for j in range(q.shape[1]):
                colour = 'white' if z[i, j] > 1.45 else INK
                ax.text(j, i, f'{q[i, j]:.3f}', ha='center', va='center', fontsize=fontsize, color=colour)
                if q[i, j] < 0.05:
                    ax.add_patch(Rectangle((j - 0.485, i - 0.47), 0.97, 0.94, fill=False, edgecolor=COL[2], linewidth=1.1))
        ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
        ax.grid(which='minor', color='white', linewidth=1.1)
        ax.tick_params(which='minor', length=0)
        return im
    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(3, 2, left=0.15, right=0.975, top=0.89, bottom=0.14, wspace=0.57, hspace=0.72, height_ratios=[1.18, 0.94, 0.92])
    a = fig.add_subplot(gs[0, 0])
    b = fig.add_subplot(gs[0, 1])
    c = fig.add_subplot(gs[1, 0])
    d = fig.add_subplot(gs[1, 1])
    e = fig.add_subplot(gs[2, 0])
    f = fig.add_subplot(gs[2, 1])
    head(a, 'A', 'Primary tests across distinct contexts', 'Directional ROAST q; BH adjustment within each contrast')
    heat(a, qmat.to_numpy(), labels, ['Culture', 'Paired', 'PRA', 'RPRA', 'Blood'], fontsize=8.8)
    a.text(0, -0.22, 'Paired and PRA share one source cohort; blood is not an airway sample.', transform=a.transAxes, fontsize=9, color=COL[5])
    head(b, 'B', 'Within-donor trajectories in nasal epithelium', 'Six paired donors; 12–18 minus 3–6 months')
    donor_names = sorted(changes.donor.unique(), key=lambda x: int(x.split()[-1]))
    y = np.arange(10)
    for j, donor in enumerate(donor_names):
        vals = changes[changes.donor == donor].set_index('gene_set').loc[programs, 'paired_logCPM_change']
        b.scatter(vals, y + (j - 2.5) * 0.075, s=19, color=COL[j], zorder=3, alpha=0.95, edgecolors='white', linewidths=0.3)
    b.scatter(summary.loc[programs, 'mean_change'], y, marker='D', s=24, c=INK, zorder=5, edgecolors='white', linewidths=0.45)
    b.axvline(0, color=COL[5], linewidth=0.8)
    b.set_yticks(y, labels)
    b.invert_yaxis()
    b.set_ylim(9.6, -0.6)
    b.set_xlabel('Mean retained-gene logCPM change (late − early)')
    polish(b)
    b.text(1, -0.29, 'Coloured points: donors; dark diamonds: donor means', transform=b.transAxes, ha='right', fontsize=9, color=COL[5])
    head(c, 'C', 'Leave-one-donor-out sensitivity', 'Four primary q < 0.05 programmes; five donors per model')
    lq = loo.pivot(index='gene_set', columns='omitted_donor', values='FDR').loc[positive, donor_names]
    clabels = []
    for g in positive:
        n = int((lq.loc[g] < 0.05).sum())
        clabels.append(label[g] + f'\n({n}/6 models at q < 0.05)')
    heat(c, lq.to_numpy(), clabels, [x.replace('patient ', 'P') for x in donor_names], fontsize=9.3)
    c.set_xlabel('Omitted donor')
    c.text(0, -0.32, 'Threshold changes also reflect the smaller sample size.', transform=c.transAxes, fontsize=9, color=COL[5])
    head(d, 'D', 'Inference depends on the test', 'Paired cohort: two procedures, separate BH adjustments')
    dq = qmat['D'].to_numpy()
    sq = summary.loc[programs, 'signflip_FDR'].to_numpy()
    for i in y:
        d.plot([dq[i], sq[i]], [i, i], color=GRID, linewidth=1.2, zorder=1)
    d.scatter(dq, y, s=22, c=COL[1], label='Directional ROAST', zorder=4)
    d.scatter(sq, y, s=24, c=COL[2], marker='D', label='Exact donor sign flips', zorder=4)
    d.set_xscale('log')
    d.set_xlim(0.009, 1)
    d.set_xticks([0.01, 0.05, 0.1, 0.5, 1], ['0.01', '0.05', '0.1', '0.5', '1'])
    d.axvline(0.05, color=COL[5], linestyle=(0, (3, 3)), linewidth=0.8)
    d.set_yticks(y, labels)
    d.set_ylim(9.6, -0.6)
    d.set_xlabel('Within-procedure q (log scale)')
    polish(d)
    d.legend(loc='lower left', bbox_to_anchor=(0, -0.46), ncol=2, frameon=False, fontsize=9, handletextpad=0.45, columnspacing=1)
    head(e, 'E', 'TGF-β sensitivity to composition', 'Cross-sectional nasal contrast: PRA versus no PRA')
    ey = np.arange(3)
    for i, row in es.iterrows():
        e.plot([row.FDR, row['FDR.Mixed']], [i, i], color=GRID, linewidth=1.2, zorder=1)
        e.scatter(row.FDR, i + 0.07, s=36, c=COL[1], zorder=3)
        e.scatter(row['FDR.Mixed'], i - 0.07, s=38, c=COL[2], marker='D', zorder=3)
        e.annotate(f"{row['FDR.Mixed']:.4f}".rstrip('0'), (row['FDR.Mixed'], i - 0.07), xytext=(0, 9), textcoords='offset points', ha='center', fontsize=9, color=COL[2])
        e.annotate(f'{row.FDR:.4f}'.rstrip('0'), (row.FDR, i + 0.07), xytext=(0, -14), textcoords='offset points', ha='center', fontsize=9, color=COL[1])
    e.set_yticks(ey, ['Unadjusted\n(n = 37)', 'Composition\n(n = 37)', 'Restricted + composition\n(n = 34)'])
    e.set_ylim(2.6, -0.65)
    e.set_xscale('log')
    e.set_xlim(0.02, 1)
    e.set_xticks([0.02, 0.05, 0.1, 0.2, 0.5, 1], ['0.02', '0.05', '0.1', '0.2', '0.5', '1'])
    e.axvline(0.05, color=COL[5], linestyle=(0, (3, 3)), linewidth=0.8)
    e.set_xlabel('Within-model q (log scale)')
    polish(e)
    e.legend(handles=[Line2D([], [], marker='o', linestyle='', color=COL[1], label='Directional'), Line2D([], [], marker='D', linestyle='', color=COL[2], label='Mixed')], loc='lower left', bbox_to_anchor=(0, -0.47), ncol=2, frameon=False, fontsize=9, handletextpad=0.5)
    head(f, 'F', 'Single-cell evidence remains donor-limited', 'Eleven donor-level nasal pseudobulks; inferred covariate linkage')
    for pos, (group, col) in enumerate([('Healthy', COL[1]), ('RPRA', COL[2])]):
        subset = donors[donors.Status == group].sort_values(['Age', 'donor'])
        offsets = np.linspace(-0.12, 0.12, len(subset))
        for j, (_, r) in enumerate(subset.iterrows()):
            f.scatter(r.Age, pos + offsets[j], s=54, color=col, marker='o' if r.male else 'D', edgecolors='white', linewidth=0.55, zorder=3)
        f.plot([subset.Age.min(), subset.Age.max()], [pos, pos], color=col, alpha=0.35, linewidth=1, zorder=1)
        f.text(66, pos, f'{int(subset.Age.min())}–{int(subset.Age.max())}', ha='right', va='center', fontsize=9, color=col)
    f.set_yticks([0, 1], ['Healthy\n(n = 6 donors)', 'RPRA\n(n = 5 donors)'])
    f.set_ylim(1.65, -0.65)
    f.set_xlim(20, 67)
    f.set_xticks([20, 30, 40, 50, 60])
    f.set_xlabel('Donor age (years)')
    polish(f)
    f.legend(handles=[Line2D([], [], marker='o', linestyle='', color=COL[5], label='Male'), Line2D([], [], marker='D', linestyle='', color=COL[5], label='Female')], loc='lower left', bbox_to_anchor=(0, -0.47), ncol=2, frameon=False, fontsize=9, handletextpad=0.5)
    f.text(0, -0.61, 'Age overlap is limited; each point is one donor.', transform=f.transAxes, fontsize=9, color=COL[5])
    save(fig, 'Figure3')

def figure_4():
    inp = TABLES / 'clinical'
    d = pd.read_csv(inp / 'clinical_lemg_merged.csv')
    props = pd.read_csv(inp / 'lemg_descriptive_proportions.csv')
    freq = pd.read_csv(inp / 'clinical_frequency_table.csv')
    tests = pd.read_csv(inp / 'exploratory_within_case_symptom_contrasts.csv')
    summary = json.loads((inp / 'laryngeal_summary.json').read_text(encoding='utf-8'))
    d['abnormal'] = d.LEMG_abnormal.astype(str).str.lower().eq('true')
    assert len(d) == 38 and d.abnormal.sum() == 29
    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(2, 3, left=0.09, right=0.975, bottom=0.14, top=0.855, wspace=0.53, hspace=0.66)
    axs = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]
    a = axs[0]
    label(a, 'A', 'Clinical cohort and assessments')
    a.axis('off')
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    a.text(0.04, 0.92, 'Published post-COVID cough series', fontsize=11.2, fontweight='bold', color=C['navy'], va='center')
    a.text(0.04, 0.8, 'Cough present ≥12 weeks\n after COVID-19 diagnosis'.replace('\n ', '\n'), fontsize=10, color=C['gray'], va='center')
    for y, n, txt in [(0.53, '38', 'Clinical records'), (0.25, '38', 'Linked LEMG records')]:
        a.text(0.04, y, n, fontsize=26, fontweight='bold', color=C['navy'], va='center')
        a.text(0.3, y, txt, fontsize=10.5, va='center', wrap=True)
    for y in [0.67, 0.39]:
        a.annotate('', xy=(0.12, y - 0.06), xytext=(0.12, y + 0.055), arrowprops={'arrowstyle': '->', 'color': C['gray'], 'lw': 1})
    a.text(0.04, 0.04, '36 normal vocal-fold movement | 2 vocal-fold palsy', fontsize=9, color=C['gray'])
    a = axs[1]
    label(a, 'B', 'Abnormal LEMG in overlapping subsets')
    short = ['All selected cases', 'Normal cord movement', 'No respiratory disorder', 'Not hospitalised', 'Hospitalised', 'Non-smokers']
    for i, row in props.iterrows():
        y = 5 - i
        c = C['navy'] if i == 0 else C['teal']
        a.errorbar(row.fraction, y, xerr=[[row.fraction - row.wilson_95_lower], [row.wilson_95_upper - row.fraction]], fmt='o', color=c, lw=1.3, capsize=2, markersize=5)
        a.text(1.03, y, f'{row.abnormal_n}/{row.denominator}', va='center', fontsize=9)
    a.set_yticks(range(6), short[::-1])
    a.set_xlim(0, 1.22)
    a.set_xticks([0, 0.25, 0.5, 0.75, 1], ['0', '25', '50', '75', '100'])
    a.set_ylim(-0.6, 5.6)
    a.set_xlabel('Abnormal LEMG (%)')
    a.text(0.5, -0.22, 'Wilson 95% CIs; subsets overlap.', transform=a.transAxes, ha='center', fontsize=9, color=C['gray'])
    a = axs[2]
    label(a, 'C', 'Movement and LEMG are distinct')
    mat = np.array([[int(((d['Vocal cord mobility'] == 'N') & d.abnormal).sum()), int(((d['Vocal cord mobility'] == 'VC palsy') & d.abnormal).sum())], [int(((d['Vocal cord mobility'] == 'N') & ~d.abnormal).sum()), int(((d['Vocal cord mobility'] == 'VC palsy') & ~d.abnormal).sum())]])
    assert mat.tolist() == [[27, 2], [9, 0]]
    cmap = LinearSegmentedColormap.from_list('teal_counts', ['#F1F5F5', C['teal']])
    a.imshow(mat, cmap=cmap, vmin=0, vmax=27, aspect='auto', extent=[-0.5, 1.5, 1.5, -0.5])
    for i in range(2):
        for j in range(2):
            a.text(j, i, str(mat[i, j]), ha='center', va='center', fontsize=25, fontweight='bold', color='white' if mat[i, j] > 15 else C['navy'])
    a.set_xticks([0, 1], ['Normal movement\n(n = 36)', 'Vocal cord palsy\n(n = 2)'])
    a.set_yticks([0, 1], ['Abnormal\nLEMG', 'Normal\nLEMG'])
    a.tick_params(length=0)
    for spine in a.spines.values():
        spine.set_visible(False)
    a.axhline(0.5, color='white', lw=3)
    a.axvline(0.5, color='white', lw=3)
    a.text(0.5, -0.23, 'Counts, not a diagnostic-accuracy estimate.', transform=a.transAxes, ha='center', fontsize=9, color=C['gray'])
    a = axs[3]
    label(a, 'D', 'Frequently recorded features')
    sel = freq[(freq.variable != 'Cough') & (freq.yes_n >= 7)].copy().sort_values(['yes_n', 'variable'], ascending=[True, False])
    labels = {'Memory loss, brain fog, attention disorders': 'Cognitive symptoms', 'Globus pharyngis': 'Globus pharyngis', 'Dysautonomic symptoms': 'Autonomic symptoms', 'Laryngeal Edema': 'Laryngeal oedema', 'Dyspnea': 'Dyspnoea', 'Arthromyalgias': 'Arthromyalgias', 'PPIs': 'PPI use', 'GERD': 'GERD'}
    y = np.arange(len(sel))
    a.barh(y, sel.yes_n, color=[C['gray'] if v in ['PPIs', 'GERD', 'Smoker'] else C['teal'] for v in sel.variable], height=0.63)
    a.set_yticks(y, [labels.get(v, v) for v in sel.variable])
    a.set_xlim(0, 28)
    a.set_xticks([0, 10, 20])
    a.set_xlabel('Patients (of 38)')
    for yy, n in zip(y, sel.yes_n):
        a.text(n + 0.45, yy, str(n), va='center', fontsize=9)
    a.text(0.5, -0.22, 'Displayed if recorded in ≥7 cases; cough: 38/38.\nAll 22 recorded fields are retained in source data.', transform=a.transAxes, ha='center', fontsize=9, color=C['gray'])
    a = axs[4]
    label(a, 'E', 'Four-muscle coverage is incomplete')
    cols = ['Right_CT', 'Left_CT', 'Right_TA', 'Left_TA']
    missing = d[cols].eq('NP').to_numpy()
    assert missing.sum() == 12 and missing.any(axis=1).sum() == 7
    a.imshow(missing.astype(int), cmap=ListedColormap([C['pale'], C['copper']]), vmin=0, vmax=1, aspect='auto', extent=[-0.5, 3.5, 37.5, -0.5])
    a.set_xlim(-0.5, 7.3)
    a.set_ylim(37.5, -0.5)
    a.set_xticks(range(4), ['R CT', 'L CT', 'R TA', 'L TA'])
    a.set_yticks([0, 9, 19, 29, 37], [str(d.iloc[i].Patient) for i in [0, 9, 19, 29, 37]])
    a.set_ylabel('Published participant ID')
    for x in [0.5, 1.5, 2.5]:
        a.plot([x, x], [-0.5, 37.5], color='white', lw=0.7)
    a.text(4.1, 6, '12 / 152', fontsize=17, fontweight='bold', color=C['copper'])
    a.text(4.1, 10, 'Muscles not\nperformed', fontsize=10, va='top')
    a.text(4.1, 20, '7 / 38', fontsize=17, fontweight='bold', color=C['navy'])
    a.text(4.1, 24, 'Patients with\npartial examinations', fontsize=10, va='top')
    a.spines['bottom'].set_visible(False)
    a.spines['left'].set_visible(False)
    a.tick_params(length=0)
    a.legend(handles=[Patch(color=C['pale'], label='Performed'), Patch(color=C['copper'], label='Not performed')], loc='upper center', bbox_to_anchor=(0.5, -0.08), ncol=2, frameon=False, handlelength=1, fontsize=9)
    a.text(0.5, -0.24, 'CT: cricothyroid; TA: thyroarytenoid.', transform=a.transAxes, ha='center', fontsize=9, color=C['gray'])
    a = axs[5]
    label(a, 'F', 'Case-only symptom comparisons')
    a.axis('off')
    a.set_xlim(0, 1)
    a.set_ylim(0, 1)
    xs = [0.02, 0.49, 0.71, 0.94]
    for x, t, ha in zip(xs, ['Symptom', 'Present', 'Absent', 'BH q'], ['left', 'center', 'center', 'center']):
        a.text(x, 0.94, t, ha=ha, fontsize=10, fontweight='bold')
    a.text(0.6, 0.88, 'Abnormal / total LEMG', ha='center', fontsize=9, color=C['gray'])
    for i, row in tests.iterrows():
        y = 0.75 - i * 0.145
        name = 'Laryngospasm' if row.symptom == 'Laryngospam' else row.symptom
        if i % 2 == 0:
            a.axhspan(y - 0.056, y + 0.056, color='#F3F5F5', zorder=0)
        a.text(0.02, y, name, va='center', fontsize=10)
        a.text(0.49, y, f'{row.yes_abnormal}/{row.yes_abnormal + row.yes_normal}', ha='center', va='center', fontsize=11, color=C['navy'])
        a.text(0.71, y, f'{row.no_abnormal}/{row.no_abnormal + row.no_normal}', ha='center', va='center', fontsize=11, color=C['navy'])
        a.text(0.94, y, f'{row.BH_q_five_tests:.3f}', ha='center', va='center', fontsize=10)
    a.text(0.02, -0.03, 'Fisher two-sided tests; BH across five contrasts.\nLaryngospasm has a zero cell (unadjusted OR = ∞).', fontsize=9, color=C['gray'], va='top')
    save(fig, 'Figure4')

def figure_5():
    navy, teal, terra, purple, gray = ('#24485C', '#2C7C87', '#C07852', '#8C7AA7', '#5F6570')
    I = TABLES / 'clinical'
    d = pd.read_csv(I / 'trial_patient_response.csv')
    t = json.loads((I / 'trial_reanalysis.json').read_text(encoding='utf-8'))
    b = pd.read_csv(I / 'between_arm_effects.csv')
    w = pd.read_csv(I / 'within_arm_changes.csv')
    fig, axs = plt.subplots(3, 2, figsize=(14, 10.5))
    fig.subplots_adjust(left=0.15, right=0.975, top=0.895, bottom=0.12, hspace=0.8, wspace=0.5)

    def panel(ax, letter, title):
        ax.set_title(title, loc='left', pad=13, fontweight='bold', color=navy)
        ax.text(-0.16, 1.13, letter, transform=ax.transAxes, fontsize=17, fontweight='bold', color='black', va='bottom')
        ax.tick_params(length=3)
        ax.set_axisbelow(True)

    def err(ax, est, lo, hi, y, c):
        ax.errorbar(est, y, xerr=[[est - lo], [hi - est]], fmt='o', color=c, ecolor=c, elinewidth=1.8, capsize=3, ms=6)
    ax = axs[0, 0]
    panel(ax, 'A', 'Day-10 response in selected postviral cough')
    for y, (arm, c) in enumerate([('Amitriptyline', teal), ('Codeine/guaifenesin', terra)]):
        sub = d[d.arm == arm]
        for val, g in sub.groupby('response_percent_day10'):
            yy = y + (np.arange(len(g)) - (len(g) - 1) / 2) * 0.055
            ax.scatter(np.repeat(val, len(g)), yy, s=40, color=c, edgecolors='white', linewidth=0.5, zorder=3)
    ax.set_yticks([0, 1], ['Amitriptyline\nn = 15', 'Codeine/guaifenesin\nn = 13'])
    ax.set_ylim(1.42, -0.42)
    ax.set_xlim(-7, 107)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel('Self-reported improvement at day 10 (%)')
    ax.grid(axis='x', alpha=0.17)
    ax.text(0, -0.35, 'Predictable chart-number allocation; no objective cough count', transform=ax.transAxes, fontsize=9, color=gray)
    ax = axs[0, 1]
    panel(ax, 'B', 'Two response definitions, one patient sample')
    for y, r in enumerate(t['thresholds']):
        lo, hi = np.array(r['risk_difference_newcombe_wilson_95_ci']) * 100
        est = r['risk_difference'] * 100
        err(ax, est, lo, hi, y, teal)
        counts = f"{r['arms']['Amitriptyline']['responders']}/15 vs {r['arms']['Codeine/guaifenesin']['responders']}/13"
        ax.text(99, y + 0.19, f'{est:.1f} [{lo:.1f}, {hi:.1f}] pp', ha='right', fontsize=9, color=gray)
        ax.text(4, y - 0.13, counts, fontsize=9)
    ax.set_yticks([0, 1], ['Response ≥50%', 'Response ≥75%'])
    ax.set_ylim(1.5, -0.55)
    ax.set_xlim(0, 103)
    ax.axvline(0, color=gray, lw=0.8)
    ax.set_xlabel('Risk difference (percentage points), 95% CI')
    ax.text(0, -0.35, 'Newcombe–Wilson intervals; thresholds are correlated', transform=ax.transAxes, fontsize=9, color=gray)
    ax = axs[1, 0]
    panel(ax, 'C', 'Postinfectious cough: between-arm LCQ contrasts')
    labels = ['Montelukast · week 2', 'Montelukast · week 4*', 'Procaterol · common\nweek 2/4 coefficient']
    for y, r in b.iterrows():
        err(ax, r['mean'], r['lower'], r['upper'], y, teal if y < 2 else purple)
    ax.set_yticks(range(3), labels)
    ax.set_ylim(2.6, -0.6)
    ax.set_xlim(-3.1, 1.1)
    ax.axvline(0, color=gray, lw=0.8, ls='--')
    ax.set_xlabel('Adjusted LCQ difference, 95% CI\nPositive values favour active treatment')
    ax = axs[1, 1]
    panel(ax, 'D', 'Within-arm recovery in the montelukast trial')
    for arm, c, off in [('Active', teal, -0.035), ('Placebo', terra, 0.035)]:
        z = w[(w.study_id == 'MONT2014') & (w.arm == arm)].sort_values('week')
        ax.errorbar(z.week + off, z['mean'], yerr=[z['mean'] - z.lower, z.upper - z['mean']], fmt='o-', color=c, capsize=3, lw=1.5, label='Montelukast' if arm == 'Active' else 'Placebo')
    ax.set_xlim(1.65, 4.35)
    ax.set_xticks([2, 4], ['Week 2', 'Week 4*'])
    ax.set_ylim(0, 7.5)
    ax.set_ylabel('LCQ improvement from baseline')
    ax.legend(frameon=False, loc='upper left', fontsize=9)
    ax.grid(axis='y', alpha=0.17)
    ax = axs[2, 0]
    panel(ax, 'E', 'Objective cough-frequency reduction in RCC/UCC')
    gefa = pd.read_csv(I / 'gefapixant_effects.csv')
    rows = [(f'{r.trial} · week {r.week}', r.reduction_percent, r.lower, r.upper) for r in gefa.itertuples()]
    for y, (label, m, l, h) in enumerate(rows):
        err(ax, m, l, h, y, teal)
        ax.text(36, y + 0.24, f'{m:.1f}% [{l:.1f}, {h:.1f}]', ha='right', fontsize=9, color=gray)
    ax.set_yticks([0, 1], [x[0] for x in rows])
    ax.set_ylim(1.6, -0.6)
    ax.set_xlim(-2, 38)
    ax.axvline(0, color=gray, lw=0.8, ls='--')
    ax.set_xlabel('Relative reduction in 24-hour cough frequency (%)\nGefapixant 45 mg twice daily versus placebo')
    ax = axs[2, 1]
    panel(ax, 'F', 'Challenge response versus clinical outcomes')
    endpoints = pd.read_csv(I / 'trial_endpoint_matrix.csv')
    cells = endpoints.pivot(index='study', columns='endpoint', values='category').loc[['SB-705498', 'XEN-D0501', 'Gefapixant', 'Amitriptyline'], ['Capsaicin challenge', 'Daily cough frequency', 'Patient-reported outcomes']].to_numpy()
    ax.imshow(cells, cmap=ListedColormap(['#EDF0F1', '#F0DED5', '#CBE0E2']), vmin=0, vmax=2, aspect='auto')
    ax.set_xticks(range(3), ['Capsaicin\nchallenge', 'Daily cough\nfrequency', 'Patient-reported\noutcomes'])
    ax.tick_params(axis='both', length=0, labelsize=9)
    ax.set_yticks(range(4), ['SB-705498\nTRPV1', 'XEN-D0501\nTRPV1', 'Gefapixant\nP2X3', 'Amitriptyline\nSelected postviral'])
    words = {0: 'Not assessed', 1: 'No benefit\ndemonstrated', 2: 'Benefit\nreported'}
    for i in range(4):
        for j in range(3):
            ax.text(j, i, 'Challenge\nmodified' if j == 0 and cells[i, j] == 2 else 'LCQ response\nCOUGH-2' if i == 2 and j == 2 else words[cells[i, j]], ha='center', va='center', fontsize=9, color=navy)
    for sp in ax.spines.values():
        sp.set_visible(False)
    save(fig, 'Figure5')

def figure_s1():
    COL = {0.3: '#2C7C87', 0.5: '#C07852'}
    GRAY = '#5F6570'
    GRID = '#E8EDF0'
    INK = '#1E2930'
    curvepath = TABLES / 'prospective/prospective_design_curves.csv'
    tablepath = TABLES / 'prospective/prospective_design_all_scenarios.csv'
    curves = pd.read_csv(curvepath)
    design = pd.read_csv(tablepath)
    curves = curves[np.isclose(curves.cough_fraction, 0.25) & np.isclose(curves.alpha, 0.05) & curves.objective.eq('detect_two_sided')]
    selected = design[design.cough_fraction_scenario.eq('cough_25pct') & design.alpha_scheme.eq('single_test_0.05') & design.delta.isin([0.3, 0.5])]
    assert len(curves) == 118
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12.4, 5.6), gridspec_kw={'width_ratios': [1.02, 1.12]})
    fig.subplots_adjust(left=0.075, right=0.98, top=0.83, bottom=0.25, wspace=0.3)

    def title(ax, letter, text, subtitle):
        ax.set_title(text, loc='left', fontsize=13, fontweight='bold', pad=29)
        ax.text(-0.12, 1.155, letter, transform=ax.transAxes, fontsize=16, fontweight='bold', va='bottom')
        ax.text(0, 1.055, subtitle, transform=ax.transAxes, fontsize=10, color=GRAY, va='bottom')
    title(ax, 'A', 'Planning for two-sided detection', 'Target cough fraction 25%; two-sided α = 0.05')
    for target, ls in [(0.8, (0, (4, 3))), (0.9, (0, (2, 3)))]:
        ax.axhline(target, color=GRAY, ls=ls, lw=0.8, zorder=0)
        ax.text(797, target + 0.012, f'{target:.0%}', ha='right', fontsize=9, color=GRAY, bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 0.4})
    landmarks = []
    for delta in [0.3, 0.5]:
        c = curves[np.isclose(curves.delta, delta)].sort_values('n_total')
        ax.plot(c.n_total, c.planned_probability, color=COL[delta], lw=2.1, label=f'δ = {delta:.1f}')
        for target in [0.8, 0.9]:
            row = selected[selected.objective.eq('detect_two_sided') & np.isclose(selected.delta, delta) & np.isclose(selected.power_target, target)].iloc[0]
            n = int(row.n_total_analyzable)
            prob = float(row.achieved_power)
            ax.scatter(n, prob, s=30, color=COL[delta], edgecolor='white', linewidth=0.6, zorder=4)
            ax.annotate(f'N = {n}', (n, prob), xytext=(0, -18), textcoords='offset points', ha='center', fontsize=9.5, color=COL[delta], bbox={'facecolor': 'white', 'edgecolor': 'none', 'pad': 0.7})
            landmarks.append({'delta': delta, 'target_power': target, 'N': n, 'exact_power': prob})
    ax.set_xlim(0, 820)
    ax.set_ylim(0, 1.015)
    ax.set_xticks([0, 200, 400, 600, 800])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0', '20', '40', '60', '80', '100'])
    ax.set_xlabel('Analysable participants, N', labelpad=8)
    ax.set_ylabel('Prospective detection power, %', labelpad=8)
    ax.grid(axis='y', color=GRID, lw=0.55, zorder=0)
    ax.legend(frameon=False, loc='lower right', title='True standardised difference', title_fontsize=9.5)
    ax.text(0, -0.24, 'Dots: exact minimum integer sample sizes at 80% and 90% power.', transform=ax.transAxes, fontsize=9, color=GRAY, va='top')
    title(bx, 'B', 'Distinguishing design objectives', '80% target power; target cough fraction 25%')
    objectives = ['detect_two_sided', 'exclude_positive_margin_two_sided_CI', 'equivalence_tost']
    labels = ['Two-sided\ndetection', 'Upper endpoint of\n95% CI below +δ', 'TOST equivalence\nwithin ±δ']
    x = np.arange(3)
    width = 0.29
    bar_values = []
    for offset, delta in [(-width / 2, 0.3), (width / 2, 0.5)]:
        rows = [selected[selected.objective.eq(obj) & np.isclose(selected.delta, delta) & np.isclose(selected.power_target, 0.8)].iloc[0] for obj in objectives]
        ns = [int(r.n_total_analyzable) for r in rows]
        bx.bar(x + offset, ns, width=width, color=COL[delta], edgecolor='none', zorder=3)
        for xx, n, obj in zip(x + offset, ns, objectives):
            bx.text(xx, n + 10, str(n), ha='center', fontsize=10.5, fontweight='bold', color=COL[delta])
            bar_values.append({'objective': obj, 'delta': delta, 'power': 0.8, 'N': n})
    bx.set_xticks(x, labels)
    bx.tick_params(axis='x', length=0, pad=9)
    bx.set_xlim(-0.55, 2.55)
    bx.set_ylim(0, 575)
    bx.set_yticks([0, 100, 200, 300, 400, 500])
    bx.set_ylabel('Required analysable participants, N', labelpad=7)
    bx.grid(axis='y', color=GRID, lw=0.55, zorder=0)
    bx.set_axisbelow(True)
    bx.text(0, -0.24, 'True difference: δ for detection; 0 for exclusion and equivalence.', transform=bx.transAxes, fontsize=9, color=GRAY, va='top')
    fig.text(0.52, 0.025, 'Illustrative future designs; thresholds are not validated clinical margins.', ha='center', fontsize=10, fontstyle='italic', color=GRAY)
    save(fig, 'FigureS1')

def main():
    font = setup()
    for draw in (figure_1, figure_2, figure_3, figure_4, figure_5, figure_s1):
        draw()
    files = sorted((p for p in OUT.iterdir() if p.suffix in {'.tiff', '.pdf', '.svg'}))
    metadata = dict(font_family='Times New Roman', font_sha256=hashlib.sha256(Path(font).read_bytes()).hexdigest(), dpi=600, files=[dict(file=p.as_posix(), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in files])
    (OUT / 'figure_manifest.json').write_text(json.dumps(metadata, indent=2) + '\n', encoding='utf-8')
    print('Six figures exported.')
if __name__ == '__main__':
    main()
