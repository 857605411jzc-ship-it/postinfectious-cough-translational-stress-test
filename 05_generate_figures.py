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
ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'output/figures'
TABLES = ROOT / 'output/tables'
C = {'navy': '#24485C', 'teal': '#2C7C87', 'copper': '#C07852', 'purple': '#8C7AA7', 'sage': '#8B9D77', 'gray': '#5F6570', 'pale': '#E8EDF0', 'ink': '#1E2930'}

def setup():
    for path in font_manager.findSystemFonts(fontext='ttf'):
        try:
            font_manager.fontManager.addfont(path)
        except RuntimeError:
            continue
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

def figure_s1():
    fig = plt.figure(figsize=(12, 8.4))
    gs = fig.add_gridspec(2, 2, left=.08, right=.97, bottom=.09, top=.94, hspace=.48, wspace=.35)
    ax = fig.add_subplot(gs[0, 0])
    label(ax, 'A', 'Public phenotype availability')
    labels = ['INCOV and MY-LC', 'Airway sources', 'Blood annotation', 'Other molecular releases']
    values = [2, 3, 1, 5]
    ax.barh(range(4), values, color=[C['teal'], C['navy'], C['copper'], '#AAB3B7'], height=.55)
    ax.set_yticks(range(4), labels); ax.invert_yaxis()
    ax.set_xlim(0, 6); ax.set_xticks(range(7)); ax.set_xlabel('Public data sources examined')
    for i, v in enumerate(values): ax.text(v+.12, i, str(v), va='center')
    ax.text(0, -.29, 'Two sources linked explicit cough status, molecular data\nand all five model covariates (Table S1).', transform=ax.transAxes, fontsize=9.5)
    ax = fig.add_subplot(gs[0, 1]); label(ax, 'B', 'Proteomic model populations')
    ax.axis('off'); ax.set_xlim(0,1); ax.set_ylim(0,1)
    for y, cohort, n, positive, negative, days, col in [( .77,'INCOV',105,26,79,'64.5 (53.3–92.0)',C['navy']),(.30,'MY-LC',97,24,73,'433.0 (404.0–463.0)',C['teal'])]:
        ax.text(0,y,cohort,fontsize=14,fontweight='bold',color=col)
        ax.text(.34,y,f'{n} participants',fontsize=14,fontweight='bold',color=col)
        ax.text(.34,y-.12,f'{positive} cough / {negative} no cough',fontsize=11)
        ax.text(.34,y-.24,f'{days} days',fontsize=10)
    ax.text(0,-.05,'Timing is median (IQR) from documented symptom onset.\nOne eligible MY-LC participant lacked primary assays.',fontsize=9.3)
    ax = fig.add_subplot(gs[1,0]); label(ax,'C','Five transcriptomic contrasts')
    names=['Paired cultures','Longitudinal nasal','Radiological groups','Nasal pseudobulk','Blood annotation']
    counts=[6,6,37,11,17]
    ax.barh(range(5),counts,color=[C['navy'],C['teal'],C['teal'],C['copper'],C['gray']],height=.55)
    ax.set_yticks(range(5),names);ax.invert_yaxis();ax.set_xlim(0,43);ax.set_xlabel('Independent donors within each contrast')
    for i,v in enumerate(counts):ax.text(v+.5,i,str(v),va='center')
    ax.text(0,-.26,'Longitudinal and radiological contrasts share GSE261777.\nThe two pseudobulk accessions represent one source study.',transform=ax.transAxes,fontsize=9.3)
    ax=fig.add_subplot(gs[1,1]);label(ax,'D','Simulation design');ax.axis('off')
    entries=[('Empirical calibration','2 cohorts · 6 correlated proteins · 5 covariates'),('Main grid','200 scenarios × 10,000 repetitions'),('Extensions','112 scenarios × 10,000 repetitions'),('Conditional design search','96 scenarios × 20,000 repetitions')]
    for i,(title,detail) in enumerate(entries):
        y=.93-i*.235
        ax.text(0,y,title,fontsize=11,fontweight='bold',color=C['navy'])
        ax.text(0,y-.085,detail,fontsize=10)
    ax.text(0,-.12,'Total: 5,040,000 target datasets',fontsize=12,fontweight='bold',color=C['teal'])
    save(fig,'Figure_S1')

def read_airway(context, filename):
    return pd.read_csv(TABLES / 'airway' / context / filename, sep='\t')


def figure_1():
    inp = TABLES / 'proteomics'
    elig = pd.read_csv('data/processed/proteomics/eligibility.csv')
    primary = pd.read_csv(inp / 'primary_panel_results.csv')
    exploratory = pd.read_csv(inp / 'exploratory_assay_results.csv')
    prec = pd.read_csv(inp / 'primary_precision.csv')
    assert len(primary) == 12 and len(exploratory) == 586
    eligible = elig[elig.eligible.astype(str).str.lower().eq('true')].copy()
    assert eligible.groupby('cohort').size().to_dict() == {'INCOV': 105, 'MYLC': 98}
    members = pd.read_csv(inp / 'primary_model_members.csv')[['cohort', 'participant_id']].drop_duplicates()
    model_population = eligible.merge(members, on=['cohort', 'participant_id'], validate='one_to_one')
    assert model_population.groupby('cohort').size().to_dict() == {'INCOV': 105, 'MYLC': 97}
    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(2, 3, left=0.072, right=0.975, bottom=0.135, top=0.855, wspace=0.38, hspace=0.66)
    axs = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(3)]
    a = axs[0]
    label(a, 'A', 'Participant flow')
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
    label(a, 'B', 'Assessment timing')
    rng = np.random.default_rng(20482)
    for i, (cohort, cough) in enumerate([('INCOV', 0), ('INCOV', 1), ('MYLC', 0), ('MYLC', 1)]):
        d = model_population[(model_population.cohort == cohort) & (model_population.cough == cough)].days.to_numpy()
        color = C['navy'] if cohort == 'INCOV' else C['teal']
        a.scatter(i + rng.uniform(-0.2, 0.2, len(d)), d, s=13, color=color, alpha=0.57, linewidth=0, zorder=2)
        q1, med, q3 = np.quantile(d, [0.25, 0.5, 0.75])
        a.vlines(i, q1, q3, color=color, lw=3, zorder=3)
        a.hlines(med, i - 0.23, i + 0.23, color='white', lw=3.4, zorder=4)
        a.hlines(med, i - 0.23, i + 0.23, color=color, lw=1.5, zorder=5)
        a.text(i, -0.13, f'n = {len(d)}', ha='center', transform=a.get_xaxis_transform(), fontsize=9)
    a.set_xticks(range(4), ['No cough', 'Cough', 'No cough', 'Cough'])
    a.set_ylim(0, 650)
    a.set_ylabel('Days since documented symptom onset')
    a.set_yticks([0, 200, 400, 600])
    a.set_xlim(-0.5, 3.5)
    a.text(0.125, -0.24, 'INCOV', transform=a.transAxes, color=C['navy'], ha='center', fontweight='bold')
    a.text(0.79, -0.24, 'MYLC', transform=a.transAxes, color=C['teal'], ha='center', fontweight='bold')
    a.grid(axis='y', color=C['pale'], lw=0.65, zorder=0)
    a = axs[2]
    label(a, 'C', 'Adjusted protein associations')
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
    a.set_xlim(min(-0.65, primary.ci_low.min() - .04), max(1.15, primary.ci_high.max() + .04))
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
    label(a, 'F', 'Interval compatibility')
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
    save(fig, 'Figure_1')

def figure_2():
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
    save(fig, 'Figure_2')

def figure_s2():
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
    save(fig, 'Figure_S2')

import matplotlib.colors as colours
import matplotlib.font_manager as fonts
TEAL, COPPER = "#176C73", "#B86D44"
COHORT_COLOURS = {"INCOV": TEAL, "MYLC": COPPER}
PROTEINS = ("IFNG", "TNF", "IL10", "CXCL9", "CXCL10", "CXCL11")
NS = (100, 169, 466, 509, 1000)
PROFILES = ((1., 1.), (.8, 1.), (1., .95), (.8, .95), (.8, .85))
PROFILE_COLOURS = ("#174D57", "#438F98", "#A8AD83", "#C18A56", "#9A4044")


def read(path):
    return pd.read_csv(path, keep_default_na=False, na_values=[""])


def tidy(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(width=.6, length=3)
    ax.grid(axis="y", color="#E4E5E4", lw=.45, zorder=0)
    ax.set_axisbelow(True)


def heading(fig, cell, letter, title, subtitle=""):
    box = cell.get_position(fig)
    fig.text(box.x0-.027, box.y1+.012, letter, fontsize=14, fontweight="bold", va="bottom")
    fig.text(box.x0, box.y1+.015, title, fontsize=10.5, fontweight="bold", va="bottom")
    if subtitle:
        fig.text(box.x0, box.y1-.004, subtitle, fontsize=8, color="#4C5253", va="top")


def figure_3():
    for path in (Path.home()/".local/share/fonts").glob("*.TTF"):
        fonts.fontManager.addfont(str(path))
    fonts.findfont("Times New Roman", fallback_to_default=False)
    plt.rcParams.update({"font.family":"Times New Roman", "font.size":8.5,
                         "axes.labelsize":8.5, "axes.titlesize":9, "xtick.labelsize":7.5,
                         "ytick.labelsize":8, "legend.fontsize":7.2, "axes.linewidth":.65,
                         "pdf.fonttype":42, "ps.fonttype":42, "svg.fonttype":"none",
                         "mathtext.fontset":"custom", "mathtext.rm":"Times New Roman",
                         "mathtext.it":"Times New Roman:italic", "mathtext.bf":"Times New Roman:bold",
                         "mathtext.cal":"Times New Roman"})
    estimates = read(ROOT/"output/ademp/simulation_estimates.csv")
    panels = read(ROOT/"output/ademp/simulation_panels.csv")
    calibration = read(ROOT/"output/ademp/empirical_calibration.csv")
    design = read(ROOT/"output/design/Table_2_conditional_design.csv")
    core = estimates[(estimates.extension=="main") & (estimates.method=="a_joint") & (estimates.scale=="native")]
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(11.7, 14.0), facecolor="white")
    outer = fig.add_gridspec(3, 2, left=.08, right=.975, bottom=.055, top=.95,
                            wspace=.22, hspace=.37, height_ratios=(.88,1.12,1.0))
    source_data = []

    cell = outer[0,0]
    heading(fig,cell,"A","Empirical variance and covariate structure")
    sub = cell.subgridspec(3,1,height_ratios=(1,1,.25),hspace=.8)
    for i, (column,label,cmap,normal) in enumerate((
        ("residual_variance_ratio","Residual variance ratio: cough / no cough", "RdBu_r", colours.TwoSlopeNorm(vcenter=0,vmin=-5,vmax=5)),
        ("nuisance_partial_R2",r"Protein variance explained by covariates, partial $R^2$", "GnBu", colours.Normalize(0,.30)))):
        ax=fig.add_subplot(sub[i]);matrix=np.array([calibration[calibration.cohort==c].set_index("protein").loc[list(PROTEINS),column] for c in ("INCOV","MYLC")])
        im=ax.imshow(np.log2(matrix) if i==0 else matrix, cmap=cmap,norm=normal,aspect="auto")
        ax.set_xticks(range(6),PROTEINS);ax.set_yticks(range(2),["INCOV","MY-LC"])
        ax.set_title(label,loc="left",pad=7)
        for r in range(2):
            for j in range(6):
                value=matrix[r,j]
                text=f"{value:.2f}" if i==0 else f"{100*value:.1f}%"
                ax.text(j,r,text,ha="center",va="center",fontsize=8,color="white" if (i==0 and abs(np.log2(value))>3) or (i==1 and value>.22) else "#202729")
        ax.tick_params(length=0,pad=5)
        for spine in ax.spines.values():spine.set_visible(False)
    ax=fig.add_subplot(sub[2]);ax.axis("off")
    r2=calibration.groupby("cohort").cough_R2_Z.first()
    ax.text(0,.5,fr"$R^2_{{C^*|Z}}$: INCOV {r2['INCOV']:.3f}; MY-LC {r2['MYLC']:.3f}",va="center")
    source_data.append(calibration.assign(panel="A"))

    cell=outer[0,1]
    heading(fig,cell,"B","Preprocessing and coefficient scale",r"$\delta=0.30$; N=466; perfect phenotype; source-interaction scenario")
    content=cell.subgridspec(2,1,height_ratios=(.16,1),hspace=0)[1]
    sub=content.subgridspec(1,2,wspace=.18)
    selected=estimates[(estimates.extension=="source_interaction")&(estimates.n==466)&(estimates.se==1)&(estimates.sp==1)&(estimates.effect=="dense03")]
    methods=("a_joint","b_readjusted","c_residual_groups","d_source")
    labels=("a  Joint model","b  Residual + Z, H","c  Residual only","d  Source; H omitted")
    for col,scale in enumerate(("native","subgroup_standardised")):
        ax=fig.add_subplot(sub[col]);tidy(ax)
        for m,method in enumerate(methods):
            for cohort,offset in (("INCOV",-.12),("MYLC",.12)):
                d=selected[(selected.method==method)&(selected.scale==scale)&(selected.cohort==cohort)]
                y=3-m+offset;colour=COHORT_COLOURS[cohort]
                values=d.mean_beta.to_numpy();truth=d.mean_scientific_target.to_numpy()
                ax.plot([values.min(),values.max()],[y,y],color=colour,lw=1.1)
                ax.scatter(values,np.full(6,y),s=14,color=colour,alpha=.7,zorder=4)
                ax.scatter([truth.mean()],[y],marker="|",s=90,color="#121819",lw=1.2,zorder=5)
        ax.set_ylim(-.5,3.75);ax.set_yticks(range(4),labels[::-1] if col==0 else [""]*4)
        ax.set_xlabel("Mean coefficient")
        ax.set_title("Reference SD" if col==0 else "Subgroup observed SD",pad=6)
        if col==0:ax.axvline(.3,color="#8E9798",lw=.7,ls="--")
        ax.set_xlim(.21,.38);ax.set_xticks([.22,.26,.30,.34,.38])
    box=cell.get_position(fig)
    fig.text(box.x0,box.y0-.035,"Teal: INCOV   Copper: MY-LC   |: scale-specific target\nPoints: six proteins; ranges describe protein heterogeneity.",fontsize=7.4,va="top")
    source_data.append(selected.assign(panel="B"))

    cell=outer[1,0]
    heading(fig,cell,"C","Finite-sample calibration of HC3 + t", "Perfect phenotype; six correlated protein tests")
    content=cell.subgridspec(2,1,height_ratios=(.12,1),hspace=0)[1]
    sub=content.subgridspec(2,2,wspace=.4,hspace=.72,height_ratios=(1,1))
    ax=fig.add_subplot(sub[0,0]);tidy(ax)
    for cohort,colour in COHORT_COLOURS.items():
        d=core[(core.cohort==cohort)&(core.se==1)&(core.sp==1)&(core.effect=="null")]
        g=d.groupby("n").coverage_95_target.agg(["mean","min","max"]).loc[list(NS)]*100
        ax.plot(NS,g['mean'],color=colour,marker="o",ms=3)
        ax.fill_between(NS,g['min'],g['max'],color=colour,alpha=.13)
    ax.axhline(95,color="#555D5F",ls="--",lw=.7);ax.set_ylim(90,100);ax.set_title("95% interval coverage")
    ax.set_ylabel("Coverage (%)");ax.set_xticks([100,466,1000]);ax.set_xlabel("N per cohort")
    ax=fig.add_subplot(sub[0,1]);tidy(ax)
    for cohort,colour in COHORT_COLOURS.items():
        for effect,style in (("null","-"),("sparse03","--")):
            d=panels[(panels.extension=="main")&(panels.method=="a_joint")&(panels.cohort==cohort)&(panels.se==1)&(panels.sp==1)&(panels.effect==effect)].sort_values("n")
            ax.errorbar(d.n,d.fdr*100,yerr=1.96*d.fdr_mcse*100,color=colour,ls=style,marker="o",ms=3,lw=1,capsize=2)
    ax.axhline(5,color="#555D5F",ls=":",lw=.7);ax.set_ylim(0,8);ax.set_title("BH false discovery rate")
    ax.set_ylabel("FDR (%)");ax.set_xticks([100,466,1000]);ax.set_xlabel("N per cohort")
    ax.text(.98,.96,"Solid: all-null\nDashed: one signal",ha="right",va="top",transform=ax.transAxes,fontsize=7)
    ax=fig.add_subplot(sub[1,:]);tidy(ax)
    for cohort,colour in COHORT_COLOURS.items():
        d=core[(core.cohort==cohort)&(core.se==1)&(core.sp==1)&(core.effect=="sparse03")&(core.protein=="IFNG")].sort_values("n")
        ax.errorbar(d.n,d.discovery_bh6*100,yerr=1.96*d.discovery_bh6_mcse*100,color=colour,lw=1.3,marker="o",ms=3,label=cohort.replace("MYLC","MY-LC"),capsize=2)
    ax.set_ylim(0,100);ax.set_xticks([100,466,1000]);ax.set_xlabel("N per cohort");ax.set_ylabel("BH discovery power (%)")
    ax.set_title(r"Sparse alternative: IFNG $\delta=0.30$; five null proteins",loc="left");ax.legend(frameon=False,loc="lower right")

    cell=outer[1,1]
    heading(fig,cell,"D","Phenotype quality and discovery power",r"$\delta=0.30$ for six proteins; mean BH discovery probability")
    content=cell.subgridspec(2,1,height_ratios=(.12,1),hspace=0)[1]
    sub=content.subgridspec(2,2,height_ratios=(1,.2),wspace=.3,hspace=.23)
    handles=[]
    for col,cohort in enumerate(("INCOV","MYLC")):
        ax=fig.add_subplot(sub[0,col]);tidy(ax)
        for (se,sp),colour in zip(PROFILES,PROFILE_COLOURS):
            d=panels[(panels.extension=="main")&(panels.method=="a_joint")&(panels.effect=="dense03")&(panels.cohort==cohort)&(panels.se==se)&(panels.sp==sp)].sort_values("n")
            line,=ax.plot(d.n,100*d.true_discovery_fraction,color=colour,lw=1.3,marker="o",ms=3,label=f"{se:.2f}/{sp:.2f}")
            if col==0:handles.append(line)
        ax.axhline(80,color="#8E9798",ls=":",lw=.7);ax.set_ylim(0,100)
        ax.set_xticks([100,466,1000]);ax.set_xlabel("N per cohort");ax.set_title(cohort.replace("MYLC","MY-LC"))
        if col==0:ax.set_ylabel("Discovery power (%)")
    ax=fig.add_subplot(sub[1,:]);ax.axis("off");ax.legend(handles=handles,loc="upper center",ncol=3,frameon=False,title="Sensitivity / specificity",title_fontsize=8)

    cell=outer[2,0]
    heading(fig,cell,"E","Latent-effect coverage and exclusion",r"True $\delta=0.50$; sensitivity 0.80; specificity 0.85")
    content=cell.subgridspec(2,1,height_ratios=(.13,1),hspace=0)[1]
    sub=content.subgridspec(1,2,wspace=.37)
    ax1,ax2=fig.add_subplot(sub[0]),fig.add_subplot(sub[1])
    for ax in (ax1,ax2):tidy(ax);ax.set_xticks([100,466,1000]);ax.set_xlabel("N per cohort")
    for cohort,colour in COHORT_COLOURS.items():
        d=core[(core.cohort==cohort)&(core.se==.8)&(core.sp==.85)&(core.effect=="dense05")]
        g=d.groupby("n").coverage_95_target.agg(["mean","min","max"]).loc[list(NS)]*100
        ax1.plot(NS,g['mean'],color=colour,marker="o",ms=3,lw=1.3,label=cohort.replace("MYLC","MY-LC"))
        ax1.fill_between(NS,g['min'],g['max'],color=colour,alpha=.13)
        for metric,style,marker in (("exclusion_030","-","o"),("equivalence_030","--","s")):
            gg=d.groupby("n")[metric].mean().loc[list(NS)]*100
            ax2.plot(NS,gg,color=colour,ls=style,marker=marker,ms=3,lw=1.3)
    ax1.set_ylim(0,100);ax1.axhline(95,color="#8E9798",ls=":",lw=.7);ax1.set_ylabel("Latent-effect coverage (%)");ax1.legend(frameon=False,loc="upper right")
    ax1.set_title("95% interval coverage")
    ax2.set_ylim(0,20);ax2.set_ylabel("Incorrect scientific conclusion (%)");ax2.set_title("Exclusion / equivalence")
    ax2.text(.03,.97,"Solid: exclude +0.30 SD\nDashed: TOST within ±0.30 SD",transform=ax2.transAxes,ha="left",va="top",fontsize=7)

    cell=outer[2,1]
    heading(fig,cell,"F","Conditional recruitment for 80% discovery power",r"$\delta=0.30$; six-protein BH correction; calibrated joint model")
    content=cell.subgridspec(2,1,height_ratios=(.13,1),hspace=0)[1]
    ax=fig.add_subplot(content);tidy(ax)
    profiles=((1.,1.),(.8,.95),(.8,.85));yticks=[]
    for ci,cohort in enumerate(("INCOV","MYLC")):
        for pi,(se,sp) in enumerate(profiles):
            row=5-ci*3-pi;yticks.append((row,f"{cohort.replace('MYLC','MY-LC')}  {se:.2f}/{sp:.2f}"))
            d=design[(design.cohort==cohort)&(design.se==se)&(design.sp==sp)&(design.effect=="dense03")&(design.objective=="discovery_bh6")]
            values=pd.to_numeric(d.estimated_N80,errors="coerce").dropna().to_numpy();colour=COHORT_COLOURS[cohort]
            lo=pd.to_numeric(d.N_below_80_MC95,errors="coerce").min()
            hi=pd.to_numeric(d.N_above_80_MC95,errors="coerce").max()
            if np.isfinite(lo) and np.isfinite(hi):
                ax.plot([lo,hi],[row,row],color="#CBD0D0",lw=6,alpha=.65,solid_capstyle="butt",zorder=1)
            if len(values):
                ax.plot([values.min(),values.max()],[row,row],color=colour,lw=1.3)
                ax.scatter(values,np.full(len(values),row),s=25,facecolor=colour,edgecolor="white",lw=.5,zorder=4)
                ax.text(values.max()*1.08,row,f"{int(values.min())}–{int(values.max())}",va="center",fontsize=8)
    ax.set_yticks([x[0] for x in yticks],[x[1] for x in yticks]);ax.set_ylim(-.6,5.7)
    ax.set_xscale("log");ax.set_xlim(100,5000);ax.set_xticks([100,250,500,1000,2500,5000],labels=["100","250","500","1000","2500","5000"])
    ax.axvline(466,color="#8E9798",ls="--",lw=.8);ax.text(466,5.45,"466",ha="center",va="bottom",fontsize=7,color="#555D5F")
    ax.set_xlabel("Evaluable N per cohort; range across six proteins")
    source_data.append(design.assign(panel="F"))
    panel_data=panels[(panels.extension=="main")&(panels.method=="a_joint")]
    pd.concat([core.assign(panel="C_D_E",record_type="protein"),panel_data.assign(panel="C_D",record_type="panel")],ignore_index=True).to_csv(OUT/"Figure_3_source_data.csv",index=False)
    pd.concat(source_data,ignore_index=True).to_csv(OUT/"Figure_3_calibration_and_design.csv",index=False)
    fig.savefig(OUT/"Figure_3.pdf",dpi=600)
    fig.savefig(OUT/"Figure_3.png",dpi=300)
    fig.savefig(OUT/"Figure_3.tiff",dpi=600,pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(OUT/"Figure_3.svg")
    plt.close(fig)



FIGURE_LEGENDS = {1: 'Protein associations and interval compatibility. (A) Participant flow into the 105 INCOV and 97 MY-LC model populations. (B) Assessment timing by recorded cough; marks show medians and interquartile ranges. Timing begins at documented symptom onset and does not establish continuous cough duration. (C) Adjusted protein differences with HC3 t-based 95% confidence intervals; outcomes were standardised within complete-model populations. Models include age, sex, BMI, log follow-up and acute hospitalisation. BH correction covers twelve comparisons. (D, E) Exploratory assays in INCOV (448) and MY-LC (138); dashed lines indicate nominal P=0.05. No assay met its within-cohort BH threshold. (F) Nominal and Bonferroni family-wise intervals containing positive reference values of 0.2, 0.5 and 0.8 SD. These references describe precision, not clinical importance.', 2: 'Transcriptomic programmes and donor-level sensitivity. (A) Directional ROAST q values for ten programmes in five contrasts. Each contrast has a separate ten-test family; outlined cells indicate q<0.05. (B) Six paired donors’ programme-mean changes between 3–6 and 12–18 months; diamonds indicate means. (C) Longitudinal donor-omission results; the four programmes retained significance in 2/6, 3/6, 2/6 and 1/6 models. (D) Primary ROAST and unweighted exact sign-flip q values; all 64 sign patterns were enumerated. (E) Directional and mixed-direction TGF-β results in the radiological comparison, before and after composition adjustment or deconvolution-error restriction. (F) Age and sex of eleven nasal pseudobulk donors. PRA denotes pulmonary radiological abnormalities; RPRA denotes respiratory PASC with radiographic abnormalities. The two GSE261777 contrasts share a source cohort. Full numerical results appear in Tables S7–S9.', 3: 'Empirically calibrated ADEMP simulations. (A) Residual variance ratios and covariate R²; phenotype R² uses the 25% calibration distribution. (B) Four preprocessing methods at N=466 and δ=0.30 SD in the source-interaction extension; H denotes the source interaction. Points represent six proteins; black marks indicate scale-specific targets. (C) HC3 coverage, BH false discovery rate and sparse-alternative discovery with perfect phenotyping. Error bars show ±1.96 Monte Carlo standard errors; shaded ranges span proteins. (D) Mean BH discovery for six +0.30-SD effects by phenotype sensitivity/specificity. (E) Latent-effect coverage and erroneous exclusion/equivalence for true +0.50-SD effects at Se/Sp=0.80/0.85. (F) Conditional per-protein N80 for six +0.30-SD effects; grey spans show evaluated pointwise Monte Carlo brackets. The dashed line marks the independent-group reference N=466. Means across proteins preserve within-repetition dependence; ranges are not confidence intervals. All diagnostic-error and interaction parameters are scenarios.'}

def main():
    setup()
    figure_1()
    figure_2()
    figure_s1()
    figure_s2()
    figure_3()
    for number, legend in FIGURE_LEGENDS.items():
        (OUT / f"Figure_{number}_legend.txt").write_text(legend + "\n", encoding="utf-8")


if __name__ == '__main__':
    main()
