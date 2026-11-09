"""Circulating-factor models and interval compatibility analyses."""

import argparse
import hashlib
import json
import platform
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import stats


RAW = Path("data/raw/proteomics")
PROCESSED = Path("data/processed/proteomics")
OUT = Path("output/tables/proteomics")
COMMIT = "e34c6de8583983074251e320acc7fd904c009739"
REPOSITORY = "aifimmunology/PASC-proteomics-talla-vasaikar-et-al"
GITHUB = f"https://raw.githubusercontent.com/{REPOSITORY}/{COMMIT}/Input%20data/"
SOURCES = (
    (
        "JimHeath_proteomics_metadata.xlsx",
        GITHUB + "JimHeath_proteomics_metadata.xlsx",
        "600aeec461cefa92f4a3b4c4624d2d31a1b3f2d12f8e9f5436ac9f9227da16b2",
    ),
    (
        "JimHeath_proteomics.xlsx",
        GITHUB + "JimHeath_proteomics.xlsx",
        "611a40a9dd2753c81975a15ffd4e6c8f8d535a492bcefca4ac2548122ca4d828",
    ),
    (
        "41586_2023_6651_MOESM4_ESM.xlsx",
        "https://media.springernature.com/full/springer-static/esm/"
        "art%3A10.1038%2Fs41586-023-06651-y/MediaObjects/"
        "41586_2023_6651_MOESM4_ESM.xlsx",
        "523e60d53347eb3d3cbb98f78c3af7bb4db168a1d9b399fdbcafbc1a06e3ee86",
    ),
)
PROTEINS = ("IFNG", "TNF", "IL10", "CXCL9", "CXCL10", "CXCL11")
MYLC_KEYS = ("IFN\u03b3", "TNF\u03b1", "IL-10", "CXCL9", "CXCL10", "CXCL11")
VARIANTS = ("rank_normal", "day56", "lung_restricted", "dyspnea_adjusted")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def obtain_sources(download):
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url, expected in SOURCES:
        path = RAW / name
        if not path.exists():
            if not download:
                raise FileNotFoundError(path)
            with urllib.request.urlopen(url, timeout=120) as response:
                body = response.read()
            require(hashlib.sha256(body).hexdigest() == expected,
                    f"Source checksum mismatch: {name}")
            path.write_bytes(body)
        require(digest(path) == expected, f"Source checksum mismatch: {name}")


def assay_label(value):
    return value.replace("\u00ce\u00b3", "\u03b3").replace("\u00ce\u00b1", "\u03b1")


def numeric(values):
    return pd.to_numeric(values, errors="coerce")


def yes_no(values):
    values = values.astype("string").str.strip().str.lower()
    require(set(values.dropna()) <= {"yes", "no"}, "Unexpected yes/no coding")
    return values.map({"yes": 1.0, "no": 0.0})


def prepare_incov():
    symptoms = pd.read_excel(RAW / SOURCES[0][0], sheet_name="S1.3 PASC data")
    clinical = pd.read_excel(RAW / SOURCES[0][0], sheet_name="S1.1 INCOV")
    proteome = pd.read_excel(RAW / SOURCES[1][0], sheet_name="S2.1 Proteomics")
    require(symptoms["Study subject ID"].is_unique, "Duplicate symptom identifier")
    require(clinical["Blood draw"].is_unique, "Duplicate clinical draw identifier")
    acute = clinical.loc[clinical["Blood draw"].str.endswith(("-T1", "-T2"))].copy()
    who = acute["Who Ordinal Scale"].astype("string").str.strip()
    category = numeric(who)
    require(set(who[category.isna()].dropna()) <= {"1 or 2"},
            "Unexpected acute WHO category")
    acute["hospitalised"] = category.ge(3).astype(float).where(category.notna())
    acute.loc[who.eq("1 or 2").fillna(False), "hospitalised"] = 0.0
    admission = acute.groupby("Study subject ID")["hospitalised"].max()
    t3 = proteome.loc[
        proteome["Healthy or INCOV"].eq("INCOV")
        & proteome["Blood Draw"].str.endswith("-T3")
    ].copy()
    require(t3["Blood Draw"].is_unique and t3["Patient Subject ID"].is_unique,
            "Duplicate third-visit participant or draw")
    t3 = t3.merge(
        clinical[["Blood draw", "Observation days since onset of symptoms",
                  "Asthma", "Chronic obstructive pulmonary disease", "Cigarette smoking"]],
        left_on="Blood Draw", right_on="Blood draw", how="left", validate="one_to_one",
    ).merge(
        symptoms[["Study subject ID", "cough", "shortness_of_breath"]],
        left_on="Patient Subject ID", right_on="Study subject ID",
        how="left", validate="one_to_one",
    )
    smoking = t3["Cigarette smoking"].astype("string").str.strip().str.lower()
    require(set(smoking.dropna()) <= {"never", "former", "current", "unknown"},
            "Unexpected smoking category")
    require(set(t3.sex.dropna()) <= {"Female", "Male"}, "Unexpected sex coding")
    d = pd.DataFrame({
        "participant_id": t3["Patient Subject ID"], "sample_id": t3["Blood Draw"],
        "cough": numeric(t3.cough), "age": numeric(t3.age), "bmi": numeric(t3.BMI),
        "days": numeric(t3["Observation days since onset of symptoms"]),
        "dyspnea": numeric(t3.shortness_of_breath), "asthma": yes_no(t3.Asthma),
        "copd": yes_no(t3["Chronic obstructive pulmonary disease"]),
        "smoking": smoking.map({"never": 0.0, "former": 1.0, "current": 1.0}),
        "sex_indicator": t3.sex.map({"Female": 1.0, "Male": 0.0}),
        "acute_severity": t3["Patient Subject ID"].map(admission),
    })
    d["no_lung_or_smoke"] = (d.asthma.eq(0) & d.copd.eq(0) & d.smoking.eq(0)).astype(int)
    columns = list(proteome.columns[6:])
    require(len(columns) == 454 and len(set(columns)) == 454, "INCOV assay universe differs")
    return pd.concat([d, t3[columns].apply(numeric)], axis=1), columns


def prepare_mylc():
    k = pd.read_excel(RAW / SOURCES[2][0])
    require(k.x0_LC_ID.is_unique, "Duplicate MY-LC participant identifier")
    fields = {
        "cohort_code": "x0_Censor_Cohort_ID", "original_exclusion": "x0_Censor_Complete",
        "cough": "x0_LC_Symptom_cough_ML_5", "age": "x0_Demographics_Age",
        "bmi": "x0_Demographics_BMI", "days": "x0_Acute_COVID19_DFSO",
        "acute_severity": "x0_Acute_COVID19_Hosp", "dyspnea": "x0_LC_Symptom_dyspnea_ML_5",
        "asthma": "x0_PMH_asthma", "copd": "x0_PMH_copd", "other_lung": "x0_PMH_otherlung",
        "oral_steroid": "x0_Censor_Oral_Steroid", "immunosuppressant": "x0_Censor_Immuno_Supress_Med",
    }
    d = pd.DataFrame({"participant_id": k.x0_LC_ID, "sample_id": k.x0_LC_ID})
    for target, source in fields.items():
        d[target] = numeric(k[source])
    sex = numeric(k.x0_Demographics_Sex)
    require(set(sex.dropna()) <= {1.0, 2.0}, "Unexpected MY-LC sex coding")
    d["sex_indicator"] = sex.map({1.0: 0.0, 2.0: 1.0})
    d["no_lung_or_smoke"] = (
        d.asthma.eq(0) & d.copd.eq(0) & d.other_lung.eq(0)
        & d.oral_steroid.eq(0) & d.immunosuppressant.eq(0)
    ).astype(int)
    columns = [c for c in k.columns if c.startswith("x1_Cytokines_")]
    require(len(columns) == 144 and len(set(columns)) == 144, "MY-LC assay universe differs")
    return pd.concat([d, k[columns].apply(numeric)], axis=1), columns


def eligibility(d, cohort):
    covariates = ["age", "bmi", "sex_indicator", "acute_severity"]
    reasons = []
    for _, row in d.iterrows():
        exclusions = []
        if cohort == "MYLC" and row.cohort_code != 3:
            exclusions.append("outside_long_covid_cohort")
        if cohort == "MYLC" and row.original_exclusion != 0:
            exclusions.append("source_exclusion")
        if row.cough not in (0, 1):
            exclusions.append("cough_not_documented_as_binary")
        if not np.isfinite(row.days) or row.days < 28:
            exclusions.append("assessment_before_day28_or_missing")
        if not np.isfinite(row[covariates].to_numpy(dtype=float)).all():
            exclusions.append("incomplete_covariates")
        reasons.append(";".join(exclusions))
    result = d[["participant_id", "sample_id", "cough", "days"]].copy()
    result.insert(0, "cohort", cohort)
    result["exclusion_reason"] = reasons
    result["eligible"] = result.exclusion_reason.eq("")
    return result


def bh(values):
    p = np.asarray(values, dtype=float)
    require(np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all(), "Invalid P value family")
    order = np.argsort(p, kind="stable")
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    result = np.empty(len(p))
    result[order] = np.minimum.accumulate(adjusted[::-1])[::-1].clip(0, 1)
    return result


def fit_model(d, column, variant="primary", omitted=None):
    if variant == "day56":
        d = d.loc[d.days.ge(56)]
    elif variant == "lung_restricted":
        d = d.loc[d.no_lung_or_smoke.eq(1)]
    elif variant == "dyspnea_adjusted":
        d = d.loc[d.dyspnea.isin([0, 1])]
    if omitted is not None:
        d = d.loc[d.participant_id.ne(omitted)]
    eligible_n = len(d)
    raw = numeric(d[column]).replace([np.inf, -np.inf], np.nan)
    nonpositive = int(raw.le(0).sum())
    observed = raw.notna()
    d = d.loc[observed]
    y = raw.loc[observed].to_numpy(dtype=float)
    n = len(d)
    n1 = int(d.cough.eq(1).sum())
    coverage = n / eligible_n if eligible_n else 0.0
    require(n >= 40 and min(n1, n - n1) >= 10 and coverage >= 0.8,
            f"Insufficient observations or assay coverage: {column}, {variant}")
    if variant == "rank_normal":
        y = stats.norm.ppf((stats.rankdata(y, method="average") - 0.5) / n)
    sd = np.std(y, ddof=1)
    require(np.isfinite(sd) and sd > 0, f"Degenerate assay variance: {column}")
    y = (y - y.mean()) / sd
    covariates = np.column_stack([
        d.age / 10, d.sex_indicator, d.bmi / 5, np.log(d.days), d.acute_severity,
    ])
    covariates -= covariates.mean(axis=0)
    x = np.column_stack([np.ones(n), d.cough, covariates])
    if variant == "dyspnea_adjusted":
        x = np.column_stack([x, d.dyspnea - d.dyspnea.mean()])
    rank = np.linalg.matrix_rank(x)
    require(rank == x.shape[1], f"Rank-deficient design: {column}, {variant}")
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    bread = np.linalg.inv(x.T @ x)
    residual = y - x @ beta
    leverage = np.einsum("ij,jk,ik->i", x, bread, x)
    require(leverage.max() < 1 - 1e-8, "HC3 undefined at unit leverage")
    # HC3 covariance; t reference with n - rank(X) degrees of freedom.
    adjusted_residual = residual / (1 - leverage)
    covariance = bread @ ((x.T * adjusted_residual**2) @ x) @ bread
    se = float(np.sqrt(covariance[1, 1]))
    require(np.isfinite(se) and se > 0, "Invalid cough-coefficient standard error")
    df = n - rank
    halfwidth = stats.t.ppf(0.975, df) * se
    z = np.delete(x, 1, axis=1)
    cough_residual = x[:, 1] - z @ np.linalg.lstsq(z, x[:, 1], rcond=None)[0]
    denominator = cough_residual @ cough_residual
    np.testing.assert_allclose(beta[1], cough_residual @ y / denominator, rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(se, np.sqrt(np.sum((cough_residual * adjusted_residual)**2))
                               / denominator, rtol=1e-8, atol=1e-9)
    return {
        "column": column, "variant": variant, "eligible_n": eligible_n,
        "n": n, "n_cough": n1, "n_no_cough": n - n1, "usable_fraction": coverage,
        "nonpositive_count": nonpositive,
        "beta": float(beta[1]), "se": se, "ci_low": float(beta[1] - halfwidth),
        "ci_high": float(beta[1] + halfwidth),
        "p": float(2 * stats.t.sf(abs(beta[1] / se), df)), "df": int(df),
        "max_leverage": float(leverage.max()), "outcome_sd_before_scaling": float(sd),
        "max_value_tie_fraction": float(pd.Series(y).value_counts().iloc[0] / n),
    }, d.participant_id.to_list()


def precision_tables(primary):
    d = primary.copy()
    d["ci_width"] = d.ci_high - d.ci_low
    d["ci_halfwidth"] = d.ci_width / 2
    d["symmetric_bound_covering_ci95"] = np.maximum(abs(d.ci_low), abs(d.ci_high))
    critical = stats.t.ppf(1 - 0.05 / (2 * len(d)), d.df)
    d["family95_low"] = d.beta - critical * d.se
    d["family95_high"] = d.beta + critical * d.se
    d["family95_width"] = d.family95_high - d.family95_low
    d["symmetric_bound_covering_family95"] = np.maximum(abs(d.family95_low), abs(d.family95_high))
    rows = []
    for r in d.itertuples():
        for label, low, high in (
            ("nominal95", r.ci_low, r.ci_high),
            ("bonferroni_family95", r.family95_low, r.family95_high),
        ):
            for bound in (0.2, 0.5, 0.8):
                rows.append({
                    "cohort": r.cohort, "protein": r.protein, "interval": label,
                    "reference_bound": bound, "ci_low": low, "ci_high": high,
                    "positive_bound_compatible": low <= bound <= high,
                    "negative_bound_compatible": low <= -bound <= high,
                    "entire_interval_inside_symmetric_bound": low > -bound and high < bound,
                    "excludes_all_positive_effects_at_or_above_bound": high < bound,
                    "excludes_all_negative_effects_at_or_below_minus_bound": low > -bound,
                })
    return d, pd.DataFrame(rows)


def demographic_tables(panels, members):
    selected = {}
    rows = []
    for cohort, d in panels.items():
        groups = members.loc[members.cohort.eq(cohort)].groupby("protein")
        identifiers = [set(g.participant_id) for _, g in groups]
        require(len(identifiers) == 6 and all(s == identifiers[0] for s in identifiers),
                "Primary assays have different analysis populations")
        used = d.loc[d.participant_id.isin(identifiers[0])].copy()
        selected[cohort] = used
        for group, g in (("All", used), ("Cough", used.loc[used.cough.eq(1)]),
                         ("No cough", used.loc[used.cough.eq(0)])):
            base = {"cohort": cohort, "group": group, "n": len(g)}
            for variable in ("age", "bmi", "days"):
                s = g[variable].dropna()
                rows.append({**base, "variable": variable, "category": "",
                             "observed_n": len(s), "missing_n": len(g) - len(s),
                             "median": s.median(), "q25": s.quantile(0.25),
                             "q75": s.quantile(0.75), "mean": s.mean(), "sd": s.std(ddof=1)})
            categories = {
                "sex_indicator": {1: "Female", 0: "Male"},
                "acute_severity": {1: "Hospitalised", 0: "Not hospitalised"},
                "smoking": {1: "Ever", 0: "Never"},
            }
            for variable, levels in categories.items():
                if variable not in g:
                    rows.append({**base, "variable": variable, "category": "Unavailable",
                                 "observed_n": 0, "missing_n": len(g)})
                    continue
                s = g[variable].dropna()
                for value, label in levels.items():
                    count = int(s.eq(value).sum())
                    rows.append({**base, "variable": variable, "category": label,
                                 "observed_n": len(s), "missing_n": len(g) - len(s),
                                 "count": count, "percent_of_group": 100 * count / len(g),
                                 "percent_of_observed": 100 * count / len(s) if len(s) else np.nan})
    combined = pd.concat(selected.values(), ignore_index=True)
    summary = {
        "n": len(combined), "cough_n": int(combined.cough.sum()),
        "female_n": int(combined.sex_indicator.sum()),
        "female_percent": float(100 * combined.sex_indicator.mean()),
        "age_median": float(combined.age.median()),
        "age_q25": float(combined.age.quantile(0.25)),
        "age_q75": float(combined.age.quantile(0.75)),
    }
    return pd.DataFrame(rows), summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Retrieve absent public source files")
    args = parser.parse_args()
    obtain_sources(args.download)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    panels, assays, ledgers = {}, {}, []
    for cohort, prepare in (("INCOV", prepare_incov), ("MYLC", prepare_mylc)):
        d, columns = prepare()
        require(d.participant_id.is_unique, f"Duplicate participant: {cohort}")
        ledger = eligibility(d, cohort)
        d.to_csv(PROCESSED / f"{cohort}_clinical_assays.csv", index=False)
        panels[cohort] = d.loc[ledger.eligible].copy()
        assays[cohort] = columns
        ledgers.append(ledger)
    mapping = []
    for protein, key in zip(PROTEINS, MYLC_KEYS):
        column = f"{protein}_INF"
        require(column in assays["INCOV"], f"Missing primary assay: {column}")
        matches = [c for c in assays["MYLC"] if assay_label(c).startswith(f"x1_Cytokines_{key}_Obs_")]
        require(len(matches) == 1, f"Ambiguous MY-LC primary assay: {protein}")
        mapping.extend([{"cohort": "INCOV", "protein": protein, "column": column},
                        {"cohort": "MYLC", "protein": protein, "column": matches[0]}])
    mapping = pd.DataFrame(mapping)
    mapping["source_scale"] = mapping.cohort.map({
        "INCOV": "Published adjusted residual Z score",
        "MYLC": "Published standardised score",
    })
    mapping.to_csv(PROCESSED / "primary_assay_map.csv", index=False)
    pd.DataFrame([
        {"cohort": cohort, "column": column}
        for cohort, columns in assays.items() for column in columns
    ]).to_csv(PROCESSED / "all_assays.csv", index=False)
    pd.concat(ledgers, ignore_index=True).to_csv(PROCESSED / "eligibility.csv", index=False)
    primary, sensitivity, omissions, members, exploratory = [], [], [], [], []
    for row in mapping.itertuples():
        d = panels[row.cohort]
        result, identifiers = fit_model(d, row.column)
        base = {"cohort": row.cohort, "protein": row.protein}
        primary.append({**base, **result})
        for identifier in identifiers:
            members.append({**base, "participant_id": identifier})
            result, _ = fit_model(d, row.column, omitted=identifier)
            omissions.append({**base, **result, "omitted_participant": identifier})
        for variant in VARIANTS:
            result, _ = fit_model(d, row.column, variant=variant)
            sensitivity.append({**base, **result})
    primary = pd.DataFrame(primary)
    primary["q_primary12"] = bh(primary.p)
    primary["q_within_cohort6"] = primary.groupby("cohort").p.transform(bh)
    sensitivity = pd.DataFrame(sensitivity)
    sensitivity["q_within_variant12"] = sensitivity.groupby("variant").p.transform(bh)
    omissions, members = pd.DataFrame(omissions), pd.DataFrame(members)
    for cohort, d in panels.items():
        excluded = set(mapping.loc[mapping.cohort.eq(cohort), "column"])
        rows = []
        for column in assays[cohort]:
            if column not in excluded:
                result, _ = fit_model(d, column)
                rows.append({"cohort": cohort, **result})
        table = pd.DataFrame(rows)
        table["q_exploratory"] = bh(table.p)
        exploratory.append(table)
    exploratory = pd.concat(exploratory, ignore_index=True)
    require((len(primary), len(sensitivity), len(omissions), len(exploratory))
            == (12, 48, 1212, 586), "Analysis population or hypothesis family differs")
    loo_summary = []
    for (cohort, protein), g in omissions.groupby(["cohort", "protein"]):
        beta = primary.loc[primary.cohort.eq(cohort) & primary.protein.eq(protein), "beta"].iloc[0]
        loo_summary.append({"cohort": cohort, "protein": protein, "n_loo": len(g),
                            "min_beta": g.beta.min(), "max_beta": g.beta.max(),
                            "same_direction_n": int(np.sign(g.beta).eq(np.sign(beta)).sum())})
    precision, compatibility = precision_tables(primary)
    demographics, pooled = demographic_tables(panels, members)
    tables = {
        "primary_panel_results": primary, "sensitivity_results": sensitivity,
        "leave_one_out_results": omissions, "loo_summary": pd.DataFrame(loo_summary),
        "primary_model_members": members, "exploratory_assay_results": exploratory,
        "primary_precision": precision, "effect_bound_compatibility": compatibility,
        "Table1_numeric_long": demographics,
    }
    for name, table in tables.items():
        table.to_csv(OUT / f"{name}.csv", index=False)
    summary = {
        "clinical_eligibility": {c: len(d) for c, d in panels.items()},
        "primary_models": len(primary), "sensitivity_models": len(sensitivity),
        "participant_omission_models": len(omissions),
        "exploratory_family_sizes": exploratory.groupby("cohort").size().to_dict(),
        "primary_minimum_q": float(primary.q_primary12.min()), "pooled_demographics": pooled,
        "source_sha256": {name: digest(RAW / name) for name, _, _ in SOURCES},
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "pandas": pd.__version__, "scipy": scipy.__version__},
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
