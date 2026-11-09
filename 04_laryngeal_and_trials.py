"""Laryngeal motor physiology and published clinical trial endpoints."""

import hashlib
import json
import math
import platform
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats

DATA = Path("data/clinical")
OUT = Path("output/tables/clinical")
NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
MUSCLES = ("Right_CT", "Left_CT", "Right_TA", "Left_TA")
SYMPTOMS = ("Globus pharyngis", "Dysphonia", "Dysphagia", "Laryngospam", "Dyspnea")
BINARY = ("GERD", "PPIs", "Smoker", "Admission", "ICU", "Taste and smell disorders",
          "Laryngospam", "Cough", "Globus pharyngis", "Dysphagia", "Dysphonia",
          "Odynophagia", "Fatigue", "Paresthesia", "Headache",
          "Memory loss, brain fog, attention disorders", "Dysautonomic symptoms",
          "Dyspnea", "Thoracic pain", "Arthromyalgias", "Laryngeal Edema", "Salivary pooling")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write_json(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                      allow_nan=False) + "\n", encoding="utf-8")


def bh(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(q, 1)
    return result


def wilson(k, n):
    require(0 <= k <= n and n > 0, "Invalid binomial denominator")
    z = stats.norm.ppf(0.975)
    p = k / n
    denominator = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denominator
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denominator
    return [max(0.0, centre - half), min(1.0, centre + half)]


def docx_tables(path):
    with zipfile.ZipFile(path) as archive:
        root = ET.fromstring(archive.read("word/document.xml"))
    return [[["".join(t.text or "" for t in cell.findall(".//w:t", NS)).strip()
              for cell in row.findall("./w:tc", NS)]
             for row in table.findall("./w:tr", NS)]
            for table in root.findall(".//w:tbl", NS)]


def laryngeal_analysis():
    raw = pd.read_excel(DATA / "pone.0283758.s004.xlsx", sheet_name="Data")
    clinical = raw.loc[pd.to_numeric(raw.Patient, errors="coerce").notna()].copy()
    clinical.columns = clinical.columns.str.strip()
    clinical["Patient"] = clinical.Patient.astype(int)
    require(clinical.Patient.is_unique and set(clinical.Patient) == set(range(1, 39)),
            "Clinical participant identifiers do not match the published series")
    rows = [r for r in docx_tables(DATA / "pone.0283758.s002.docx")[0]
            if r and r[0].isdigit()]
    require(len(rows) == 38 and all(len(r) == 6 for r in rows), "Incomplete LEMG source")
    lemg = pd.DataFrame(rows, columns=("Patient", *MUSCLES, "LEMG_description"))
    lemg["Patient"] = lemg.Patient.astype(int)
    lemg["LEMG_abnormal"] = lemg.LEMG_description.str.strip().ne("Normal.")
    require(lemg.Patient.is_unique, "Duplicate LEMG participant")
    abnormal = [r for r in docx_tables(DATA / "pone.0283758.s003.docx")[0]
                if r and r[0].isdigit()]
    require({int(r[0]) for r in abnormal} == set(lemg.loc[lemg.LEMG_abnormal, "Patient"]),
            "Published LEMG classifications disagree between source tables")
    merged = clinical.merge(lemg, on="Patient", validate="one_to_one")
    for row in abnormal:
        person = merged.loc[merged.Patient.eq(int(row[0]))].iloc[0]
        require(person.Sex == row[1] and float(person.Age) == float(row[2]),
                "Published demographics disagree between source tables")
    for column in BINARY:
        values = merged[column].astype("string").str.strip().str.lower()
        require(values.notna().all() and set(values) <= {"yes", "no"},
                f"Invalid binary clinical field: {column}")
        merged[column] = values.map({"yes": "Yes", "no": "No"})
    require(merged.Cough.eq("Yes").all(), "Unexpected source cough phenotype")
    merged.to_csv(OUT / "clinical_lemg_merged.csv", index=False)
    lemg.to_csv(OUT / "lemg_source_extracted.csv", index=False)
    groups = {
        "all selected post-COVID cough cases": np.ones(len(merged), dtype=bool),
        "normal vocal cord mobility": merged["Vocal cord mobility"].eq("N"),
        "no recorded respiratory disorder": merged["Respiratory disorders"].eq("No"),
        "not hospitalised acutely": merged.Admission.eq("No"),
        "hospitalised acutely": merged.Admission.eq("Yes"),
        "non-smokers": merged.Smoker.eq("No"),
    }
    proportions = []
    for group, selected in groups.items():
        n = int(selected.sum())
        k = int(merged.loc[selected, "LEMG_abnormal"].sum())
        low, high = wilson(k, n)
        ci = stats.binomtest(k, n).proportion_ci(method="wilson")
        require(np.allclose([low, high], [ci.low, ci.high], rtol=0, atol=1e-12),
                "Wilson interval verification failed")
        proportions.append(dict(group=group, abnormal_n=k, denominator=n, fraction=k/n,
                                wilson_95_lower=low, wilson_95_upper=high))
    pd.DataFrame(proportions).to_csv(OUT / "lemg_descriptive_proportions.csv", index=False)
    comparisons = []
    for symptom in SYMPTOMS:
        exposed = merged[symptom].eq("Yes")
        outcome = merged.LEMG_abnormal
        a, b = int((exposed & outcome).sum()), int((exposed & ~outcome).sum())
        c, d = int((~exposed & outcome).sum()), int((~exposed & ~outcome).sum())
        odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        comparisons.append(dict(symptom=symptom, yes_abnormal=a, yes_normal=b,
                                no_abnormal=c, no_normal=d, unadjusted_odds_ratio=odds,
                                fisher_two_sided_p=p))
    comparisons = pd.DataFrame(comparisons)
    comparisons["BH_q_five_tests"] = bh(comparisons.fisher_two_sided_p)
    comparisons.to_csv(OUT / "exploratory_within_case_symptom_contrasts.csv", index=False)
    pd.DataFrame([dict(variable=c, yes_n=int(merged[c].eq("Yes").sum()),
                       no_n=int(merged[c].eq("No").sum()), missing_n=0)
                  for c in BINARY]).to_csv(OUT / "clinical_frequency_table.csv", index=False)
    incomplete = lemg[list(MUSCLES)].eq("NP")
    summary = dict(n=len(merged), LEMG_abnormal=int(merged.LEMG_abnormal.sum()),
                   normal_LEMG=int((~merged.LEMG_abnormal).sum()),
                   age_mean=float(merged.Age.mean()), age_sd=float(merged.Age.std(ddof=1)),
                   sex_counts=merged.Sex.value_counts().to_dict(),
                   normal_mobility_n=int(merged["Vocal cord mobility"].eq("N").sum()),
                   muscles_not_performed=int(incomplete.sum().sum()),
                   participants_with_incomplete_four_muscle_exam=int(incomplete.any(axis=1).sum()),
                   normal_LEMG_with_incomplete_four_muscle_exam=int((~lemg.LEMG_abnormal & incomplete.any(axis=1)).sum()))
    require(summary["LEMG_abnormal"] == 29 and summary["normal_LEMG"] == 9,
            "LEMG classification count differs from the source series")
    write_json("laryngeal_summary.json", summary)


def treatment_contrast(rows, threshold):
    counts = {}
    for arm in ("Amitriptyline", "Codeine/guaifenesin"):
        values = rows.loc[rows.arm.eq(arm), "response_percent_day10"]
        n = len(values)
        k = int(values.ge(threshold).sum())
        counts[arm] = dict(responders=k, n=n, nonresponders=n-k,
                           risk=k/n, wilson_95_ci=wilson(k, n))
    active, control = counts.values()
    a, b = active["responders"], active["nonresponders"]
    c, d = control["responders"], control["nonresponders"]
    odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    total, nactive, responders = a+b+c+d, a+b, a+c
    probability = lambda k: (math.comb(responders, k) * math.comb(total-responders, nactive-k)
                              / math.comb(total, nactive))
    support = range(max(0, nactive-(total-responders)), min(nactive, responders)+1)
    enumerated = sum(probability(k) for k in support if probability(k) <= probability(a)*(1+1e-12))
    require(abs(p-enumerated) < 1e-12, "Fisher probability verification failed")
    p1, p0 = active["risk"], control["risk"]
    l1, u1 = active["wilson_95_ci"]
    l0, u0 = control["wilson_95_ci"]
    difference = p1-p0
    # Newcombe's uncorrected hybrid-score interval combines the marginal Wilson limits.
    interval = [max(-1.0, difference-math.sqrt((p1-l1)**2+(u0-p0)**2)),
                min(1.0, difference+math.sqrt((u1-p1)**2+(p0-l0)**2))]
    ratio = ratio_interval = logodds = logodds_se = wald_p = None
    if min(a, b, c, d) > 0:
        ratio = p1/p0
        se = math.sqrt(1/a-1/(a+b)+1/c-1/(c+d))
        z = stats.norm.ppf(.975)
        ratio_interval = [math.exp(math.log(ratio)-z*se), math.exp(math.log(ratio)+z*se)]
        logodds = math.log(a*d/(b*c))
        logodds_se = math.sqrt(1/a+1/b+1/c+1/d)
        wald_p = float(2*stats.norm.sf(abs(logodds/logodds_se)))
    return dict(threshold_percent=threshold, arms=counts,
                table_active_control_rows_responder_nonresponder_cols=[[a,b],[c,d]],
                fisher_exact_two_sided_p=float(p), enumerated_fisher_p=enumerated,
                sample_odds_ratio=float(odds) if math.isfinite(odds) else None,
                odds_ratio_unbounded=not math.isfinite(odds),
                risk_difference=difference, risk_difference_newcombe_wilson_95_ci=interval,
                risk_ratio=ratio, risk_ratio_log_wald_95_ci=ratio_interval,
                log_odds_ratio=logodds, log_odds_ratio_standard_error=logodds_se,
                log_odds_ratio_wald_two_sided_p=wald_p)


def trial_analysis():
    rows = pd.read_csv(DATA / "trial_patient_response.csv")
    require(rows.published_patient_number.tolist() == list(range(1, 29)),
            "Incomplete published participant-response table")
    require(rows.response_percent_day10.isin([0,25,50,75,100]).all(), "Invalid response category")
    require(rows.arm.value_counts().to_dict() == {"Amitriptyline":15, "Codeine/guaifenesin":13},
            "Unexpected trial arm denominators")
    results = [treatment_contrast(rows, t) for t in (50,75)]
    require(results[0]["table_active_control_rows_responder_nonresponder_cols"] == [[13,2],[1,12]]
            and results[1]["table_active_control_rows_responder_nonresponder_cols"] == [[11,4],[0,13]],
            "Response thresholds disagree with the published table")
    rows.to_csv(OUT / "trial_patient_response.csv", index=False)
    write_json("trial_reanalysis.json", dict(n=len(rows), thresholds=results,
               source_doi="10.1097/01.mlg.0000244377.60334.e3", source_table="Table III"))
    between = pd.read_csv(DATA / "between_arm_effects.csv")
    within = pd.read_csv(DATA / "within_arm_changes.csv")
    procaterol = within.study_id.eq("PROC2021")
    n = within.loc[procaterol, "n_primary_analyzed"]
    margin = stats.t.ppf(.975, n-1)*within.loc[procaterol, "sd"]/np.sqrt(n)
    within.loc[procaterol, "lower"] = within.loc[procaterol, "mean"]-margin
    within.loc[procaterol, "upper"] = within.loc[procaterol, "mean"]+margin
    between.to_csv(OUT / "between_arm_effects.csv", index=False)
    within.to_csv(OUT / "within_arm_changes.csv", index=False)
    for name in ("gefapixant_effects.csv", "trial_endpoint_matrix.csv"):
        pd.read_csv(DATA / name).to_csv(OUT / name, index=False)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv("data/input_manifest.tsv", sep="\t")
    for row in manifest.loc[manifest.file.str.startswith("data/clinical/")].itertuples():
        path = Path(row.file)
        require(hashlib.sha256(path.read_bytes()).hexdigest() == row.sha256,
                f"Clinical source checksum mismatch: {path}")
    laryngeal_analysis()
    trial_analysis()
    write_json("environment.json", dict(python=platform.python_version(), numpy=np.__version__,
                                        pandas=pd.__version__, scipy=scipy.__version__))
    print("Clinical analyses complete.")


if __name__ == "__main__":
    main()
