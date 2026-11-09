"""Prospective detection, exclusion, equivalence and precision calculations."""

import argparse
import json
import math
import platform
from fractions import Fraction
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import integrate, special, stats


OUT = Path("output/tables/prospective")
FRACTIONS = (
    ("observed_50_of_202", Fraction(50, 202)),
    ("cough_25pct", Fraction(1, 4)),
    ("balanced_50pct", Fraction(1, 2)),
)
ALPHAS = (("single_test_0.05", 0.05), ("bonferroni_12", 0.05 / 12))
OBJECTIVES = (
    "detect_two_sided", "exclude_positive_margin",
    "exclude_positive_margin_two_sided_CI", "equivalence_tost",
)


def group_sizes(n, fraction):
    n1 = (fraction.numerator * n + fraction.denominator - 1) // fraction.denominator
    return n1, n - n1


def recruitment(n1, n0):
    return math.ceil(Fraction(n1, 1) / Fraction(4, 5)), math.ceil(Fraction(n0, 1) / Fraction(4, 5))


@lru_cache(maxsize=None)
def tost_power(n1, n0, margin, alpha, true_difference=0.0):
    if min(n1, n0) < 2 or margin <= 0 or not 0 < alpha < 0.5:
        raise ValueError("Invalid equivalence design")
    df = n1 + n0 - 2
    s = math.sqrt(1 / n1 + 1 / n0)
    critical = stats.t.ppf(1 - alpha, df)
    xmax = df * (margin / (critical * s))**2
    low = stats.chi2.ppf(1e-12, df)
    high = min(xmax, stats.chi2.ppf(1 - 1e-12, df))
    if high <= low:
        return 0.0

    # Shieh (2016): independent normal mean difference and chi-square pooled variance.
    def integrand(x):
        width = margin / s - critical * math.sqrt(x / df)
        probability = (special.ndtr(width - true_difference / s)
                       - special.ndtr(-width - true_difference / s))
        return max(0.0, probability) * stats.chi2.pdf(x, df)

    # The omitted chi-square tails contribute at most 2e-12 to power.
    result, error = integrate.quad(integrand, low, high, epsabs=2e-10, epsrel=2e-10, limit=150)
    if error > 1e-8:
        raise ArithmeticError("Equivalence integration tolerance exceeded")
    return float(result)


@lru_cache(maxsize=None)
def design_power(n, fraction, delta, alpha, objective):
    n1, n0 = group_sizes(n, fraction)
    if min(n1, n0) < 2:
        return 0.0
    df = n - 2
    noncentrality = delta / math.sqrt(1 / n1 + 1 / n0)
    if objective == "detect_two_sided":
        critical = stats.t.ppf(1 - alpha / 2, df)
        return float(stats.nct.sf(critical, df, noncentrality)
                     + stats.nct.cdf(-critical, df, noncentrality))
    if objective in ("exclude_positive_margin", "exclude_positive_margin_two_sided_CI"):
        component_alpha = alpha if objective == "exclude_positive_margin" else alpha / 2
        critical = stats.t.ppf(1 - component_alpha, df)
        return float(stats.nct.cdf(-critical, df, -noncentrality))
    if objective == "equivalence_tost":
        return tost_power(n1, n0, delta, alpha)
    raise ValueError(objective)


def minimum_total(probability, target):
    low, high = 3, 8
    while probability(high) < target:
        high *= 2
        if high > 131072:
            raise ArithmeticError("Sample-size search exceeded its domain")
    while high - low > 1:
        middle = (high + low) // 2
        if probability(middle) >= target:
            high = middle
        else:
            low = middle
    achieved, previous = probability(high), probability(high - 1)
    if not previous < target <= achieved:
        raise ArithmeticError("Integer sample-size boundary is inconsistent")
    return high, achieved, previous


def width_assurance(n, fraction, halfwidth, alpha):
    n1, n0 = group_sizes(n, fraction)
    if min(n1, n0) < 2:
        return 0.0
    df = n - 2
    critical = stats.t.ppf(1 - alpha / 2, df)
    s = math.sqrt(1 / n1 + 1 / n0)
    return float(stats.chi2.cdf(df * (halfwidth / (critical * s))**2, df))


def prospective_tables():
    rows, precision, curves, benchmarks = [], [], [], []
    for (label, fraction), (scheme, alpha), delta, target, objective in product(
        FRACTIONS, ALPHAS, (0.2, 0.3, 0.5, 0.8), (0.8, 0.9), OBJECTIVES
    ):
        probability = lambda n: design_power(n, fraction, delta, alpha, objective)
        n, achieved, previous = minimum_total(probability, target)
        n1, n0 = group_sizes(n, fraction)
        r1, r0 = recruitment(n1, n0)
        rows.append({
            "objective": objective, "delta": delta,
            "assumed_true_standardized_difference": delta if objective == "detect_two_sided" else 0.0,
            "cough_fraction_scenario": label, "target_cough_fraction": float(fraction),
            "alpha_scheme": scheme, "alpha_per_hypothesis": alpha, "power_target": target,
            "n_total_analyzable": n, "n_cough_analyzable": n1, "n_no_cough_analyzable": n0,
            "actual_cough_fraction": n1 / n, "achieved_power": achieved,
            "power_at_N_minus_1": previous, "recruit_cough_20pct_loss": r1,
            "recruit_no_cough_20pct_loss": r0, "n_recruit_20pct_loss": r1 + r0,
        })
    for (label, fraction), (scheme, alpha), halfwidth, target in product(
        FRACTIONS, ALPHAS, (0.2, 0.3, 0.5), (0.8, 0.9)
    ):
        probability = lambda n: width_assurance(n, fraction, halfwidth, alpha)
        n, achieved, previous = minimum_total(probability, target)
        n1, n0 = group_sizes(n, fraction)
        r1, r0 = recruitment(n1, n0)
        precision.append({
            "cough_fraction_scenario": label, "target_cough_fraction": float(fraction),
            "alpha_scheme": scheme, "two_sided_interval_confidence": 1 - alpha,
            "halfwidth_target_in_planning_SD": halfwidth, "width_assurance_target": target,
            "n_total_analyzable": n, "n_cough_analyzable": n1, "n_no_cough_analyzable": n0,
            "achieved_width_assurance": achieved, "assurance_at_N_minus_1": previous,
            "n_recruit_20pct_loss": r1 + r0,
        })
    for fraction, delta, (_, alpha), n, objective in product(
        (Fraction(1, 4), Fraction(1, 2)), (0.3, 0.5), ALPHAS, range(40, 1201, 20), OBJECTIVES
    ):
        curves.append({"n_total": n, "cough_fraction": float(fraction), "delta": delta,
                       "alpha": alpha, "objective": objective,
                       "planned_probability": design_power(n, fraction, delta, alpha, objective)})
    for (scheme, alpha), delta in product(ALPHAS, (0.3, 0.5)):
        n, n1, n0 = 100, 25, 75
        value = design_power(n, Fraction(1, 4), delta, alpha, "detect_two_sided")
        benchmarks.append({
            "total_n": n, "cough_n": n1, "no_cough_n": n0,
            "assumed_true_standardized_difference": delta, "two_sided_alpha": alpha,
            "exact_power": value, "exact_power_percent": 100 * value, "residual_df": n - 2,
            "noncentrality": delta / math.sqrt(1 / n1 + 1 / n0),
            "critical_t": stats.t.ppf(1 - alpha / 2, n - 2), "alpha_scope": scheme,
        })
    table = pd.DataFrame(rows)
    selected = table.loc[
        table.cough_fraction_scenario.eq("cough_25pct") & table.delta.isin([0.3, 0.5])
        & table.alpha_scheme.eq("single_test_0.05") & table.power_target.eq(0.8)
        & table.objective.ne("exclude_positive_margin")
    ]
    return {
        "prospective_design_all_scenarios": table,
        "Table_S13_selected_scenarios": selected,
        "prospective_CI_width_assurance": pd.DataFrame(precision),
        "prospective_design_curves": pd.DataFrame(curves),
        "N100_detection_power": pd.DataFrame(benchmarks),
    }


def simulation_checks():
    rng = np.random.default_rng(20261028)
    rows = []
    for objective in OBJECTIVES:
        n, fraction, delta, alpha, replicates = 200, Fraction(1, 4), 0.5, 0.05, 500000
        n1, n0 = group_sizes(n, fraction)
        df = n - 2
        s = math.sqrt(1 / n1 + 1 / n0)
        true_difference = delta if objective == "detect_two_sided" else 0.0
        difference = rng.normal(true_difference, s, replicates)
        pooled_sd = np.sqrt(rng.chisquare(df, replicates) / df)
        if objective == "detect_two_sided":
            hit = abs(difference / (pooled_sd * s)) > stats.t.ppf(1 - alpha / 2, df)
        elif objective in ("exclude_positive_margin", "exclude_positive_margin_two_sided_CI"):
            component_alpha = alpha if objective == "exclude_positive_margin" else alpha / 2
            hit = difference + stats.t.ppf(1 - component_alpha, df) * pooled_sd * s < delta
        else:
            halfwidth = stats.t.ppf(1 - alpha, df) * pooled_sd * s
            hit = (difference - halfwidth > -delta) & (difference + halfwidth < delta)
        exact = design_power(n, fraction, delta, alpha, objective)
        estimated = float(hit.mean())
        mcse = math.sqrt(estimated * (1 - estimated) / replicates)
        if abs(estimated - exact) >= 5 * mcse:
            raise ArithmeticError(f"Simulation discrepancy: {objective}")
        rows.append({"objective": objective, "replicates": replicates, "seed": 20261028,
                     "exact_probability": exact, "simulated_probability": estimated,
                     "monte_carlo_se": mcse})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate", action="store_true", help="Compare exact probabilities with simulation")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    benchmark = tost_power(49, 207, 5.92 / 9.78, 0.05, 2.2 / 9.78)
    if abs(benchmark - 0.7711) >= 0.00006:
        raise ArithmeticError("Shieh (2016) numerical benchmark differs")
    tables = prospective_tables()
    if args.validate:
        tables["simulation_checks"] = simulation_checks()
    for name, table in tables.items():
        table.to_csv(OUT / f"{name}.csv", index=False)
    summary = {
        "design_scenarios": len(tables["prospective_design_all_scenarios"]),
        "precision_scenarios": len(tables["prospective_CI_width_assurance"]),
        "shieh_2016_power": benchmark, "minimum_integer_boundaries_checked": True,
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "pandas": pd.__version__, "scipy": scipy.__version__},
    }
    (OUT / "design_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
