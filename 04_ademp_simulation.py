from __future__ import annotations

import argparse
import contextlib
import hashlib
import itertools
import json
import math
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_variable] = "1"

import numpy as np
import pandas as pd
import scipy
from scipy import integrate, linalg, stats
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent
PROTEINS = ("IFNG", "TNF", "IL10", "CXCL9", "CXCL10", "CXCL11")
COHORTS = ("INCOV", "MYLC")
METHODS = ("oracle", "a_joint", "b_readjusted", "c_residual_groups", "d_source")
ERROR_PROFILES = ((1., 1.), (.8, 1.), (1., .95), (.8, .95), (.8, .85))
NS = (100, 169, 466, 509, 1000)
EFFECTS = {"null": [0.] * 6, "sparse03": [.3, 0., 0., 0., 0., 0.],
           "dense03": [.3] * 6, "dense05": [.5] * 6,
           "negative03": [-.3] * 6, "negative05": [-.5] * 6}
SEED = 2026111306
BATCH = 64
CALIBRATION = {}


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".part")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def stream(key, block=0):
    words = np.frombuffer(hashlib.sha256(key.encode()).digest()[:16], dtype="<u4")
    return np.random.default_rng(np.random.SeedSequence([SEED, *map(int, words), int(block)]))


def weighted_cov(x, w):
    mean = w @ x
    xc = x - mean
    return mean, (xc.T * w) @ xc


def prepare(source, destination):
    source = Path(source)
    assays = pd.read_csv(source / "data/processed/proteomics/primary_assay_map.csv")
    members = pd.read_csv(source / "output/tables/proteomics/primary_model_members.csv")
    references = pd.read_csv(source / "output/tables/proteomics/primary_panel_results.csv")
    result = {"proteins": PROTEINS, "covariates": ["age", "sex", "bmi", "log_days", "acute_hospitalisation"],
              "cohorts": {}, "reference_unit": "within-group SD at zero target effect and true cough fraction 0.25"}
    for cohort in COHORTS:
        input_path = source / f"data/processed/proteomics/{cohort}_clinical_assays.csv"
        d = pd.read_csv(input_path)
        selected = set.intersection(*[set(members.loc[(members.cohort == cohort) & (members.protein == p), "participant_id"]) for p in PROTEINS])
        d = d[d.participant_id.isin(selected)].sort_values("participant_id")
        c = d.cough.to_numpy(float)
        z = d[["age", "sex_indicator", "bmi", "days", "acute_severity"]].to_numpy(float)
        z[:, 3] = np.log(z[:, 3])
        z = (z - z.mean(0)) / z.std(0, ddof=1)
        columns = [assays.loc[(assays.cohort == cohort) & (assays.protein == p), "column"].item() for p in PROTEINS]
        y = d[columns].to_numpy(float)
        y = (y - y.mean(0)) / y.std(0, ddof=1)
        x = np.column_stack((np.ones(len(c)), c, z))
        bread = np.linalg.inv(x.T @ x)
        coefficients = bread @ x.T @ y
        residual = y - x @ coefficients
        leverage = np.sum((x @ bread) * x, axis=1)
        hc2 = residual / np.sqrt(1 - leverage)[:, None]
        se = np.sqrt(np.sum((residual / (1 - leverage)[:, None] * (x @ bread)[:, 1, None]) ** 2, axis=0))
        ref = references[references.cohort == cohort].set_index("protein").loc[list(PROTEINS)]
        if not (np.allclose(coefficients[1], ref.beta, atol=1e-12) and np.allclose(se, ref.se, atol=1e-12)):
            raise ValueError(f"Primary model reproduction failed: {cohort}")
        short = x[:, :2]
        e0 = y - short @ np.linalg.lstsq(short, y, rcond=None)[0]
        partial_r2 = 1 - np.sum(residual ** 2, 0) / np.sum(e0 ** 2, 0)
        w = np.where(c == 1, .25 / sum(c == 1), .75 / sum(c == 0))
        hraw = z[:, 0] * z[:, 1]
        k = np.column_stack((np.ones(len(c)), z))
        hcoef = np.linalg.solve((k.T * w) @ k, (k.T * w) @ hraw)
        hr = hraw - k @ hcoef
        hsd = math.sqrt(np.sum(w * hr ** 2))
        item = {"n": len(c), "n_cough": int(c.sum()), "source_sha256": sha(input_path),
                "gamma_direction": coefficients[2:].tolist(), "nuisance_partial_r2": partial_r2.tolist(),
                "h_coefficients": hcoef.tolist(), "h_sd": hsd,
                "primary_beta": coefficients[1].tolist(), "primary_hc3_se": se.tolist(),
                "groups": {}}
        for g in (0, 1):
            zg, eg = z[c == g], hc2[c == g]
            item["groups"][str(g)] = {"z": zg.tolist(), "residuals": (eg - eg.mean(0)).tolist(),
                                      "covariance": (eg.T @ eg / len(eg)).tolist()}
        result["cohorts"][cohort] = item
    write_json(destination, result)


def initialise(path):
    global CALIBRATION
    CALIBRATION = json.loads(Path(path).read_text())
    threadpool_limits(limits=1)


def scenario(cohort, n, se, sp, effect, extension="main", repetitions=10000):
    return {"cohort": cohort, "n": int(n), "se": float(se), "sp": float(sp), "effect": effect,
            "extension": extension, "pi": .10 if extension == "prevalence10" else .25,
            "repetitions": int(repetitions)}


def key(s, paired=False):
    fields = [s["extension"], f"n{s['n']}", f"se{round(100*s['se'])}", f"sp{round(100*s['sp'])}", s["effect"]]
    return "_".join(([] if paired else [s["cohort"]]) + fields)


def grid(repetitions=10000):
    cells = [scenario(c, n, se, sp, e, repetitions=repetitions)
             for c, n, (se, sp), e in itertools.product(COHORTS, NS, ERROR_PROFILES,
                                                       ("null", "sparse03", "dense03", "dense05"))]
    for ext in ("homoscedastic", "r2zero", "r2high", "empirical_errors", "source_interaction", "prevalence10"):
        cells.extend(scenario(c, n, se, sp, e, ext, repetitions)
                     for c, n, (se, sp), e in itertools.product(COHORTS, (100, 466),
                                                               ((1., 1.), (.8, .95)), ("null", "dense03")))
    cells.extend(scenario(c, n, se, sp, e, "negative_boundary", repetitions)
                 for c, n, (se, sp), e in itertools.product(COHORTS, (100, 466),
                                                           ((1., 1.), (.8, .95)), ("negative03", "negative05")))
    assert len(cells) == len({key(s) for s in cells}) == 312
    return cells


def model_environment(s):
    item = CALIBRATION["cohorts"][s["cohort"]]
    zs = [np.asarray(item["groups"][str(g)]["z"]) for g in (0, 1)]
    covs = [np.asarray(item["groups"][str(g)]["covariance"]) for g in (0, 1)]
    gamma = np.asarray(item["gamma_direction"])
    q = np.asarray(item["nuisance_partial_r2"])
    if s["extension"] == "r2zero":
        q = np.zeros(6)
    if s["extension"] == "r2high":
        q = np.full(6, .30)
    within = sum(weight * np.cov(z, rowvar=False, bias=True) for weight, z in zip((.75, .25), zs))
    fvar = np.sum(gamma * (within @ gamma), axis=0)
    gamma = gamma * np.sqrt(q / fvar)
    pooled = .75 * covs[0] + .25 * covs[1]
    multiplier = np.sqrt((1 - q) / np.diag(pooled))
    covs = [v * multiplier[:, None] * multiplier[None, :] for v in covs]
    if s["extension"] == "homoscedastic":
        common = .75 * covs[0] + .25 * covs[1]
        covs = [common, common]
    factors = [np.linalg.cholesky(v) for v in covs]
    empirical = []
    for g in (0, 1):
        e = np.asarray(item["groups"][str(g)]["residuals"])
        e = e - e.mean(0)
        white = np.linalg.solve(np.linalg.cholesky(e.T @ e / len(e)), e.T).T
        empirical.append(white @ factors[g].T)
    hcoef, hsd = np.asarray(item["h_coefficients"]), item["h_sd"]
    hs = [(z[:, 0] * z[:, 1] - np.column_stack((np.ones(len(z)), z)) @ hcoef) / hsd for z in zs]
    return {"z": zs, "covariance": covs, "factor": factors, "empirical": empirical,
            "gamma": gamma, "hcoef": hcoef, "hsd": hsd, "h": hs,
            "interaction": s["extension"] == "source_interaction", "delta": np.asarray(EFFECTS[s["effect"]])}


def population(env, pi, se, sp, include_h=False):
    z = np.vstack([np.column_stack((x, h)) if include_h else x for x, h in zip(env["z"], env["h"])])
    c = np.r_[np.zeros(len(env["z"][0])), np.ones(len(env["z"][1]))]
    weights = np.where(c == 1, pi / sum(c == 1), (1 - pi) / sum(c == 0))
    mean_z, variance_z = weighted_cov(z, weights)
    cov_cz = (weights * (c - pi)) @ z
    variance_c = pi * (1 - pi)
    var_c_res = variance_c - cov_cz @ np.linalg.solve(variance_z, cov_cz)
    j = se + sp - 1
    pw = se * pi + (1 - sp) * (1 - pi)
    var_w = pw * (1 - pw)
    v = pi * se * (1 - se) + (1 - pi) * sp * (1 - sp)
    var_w_res = j*j*var_c_res + v
    theta_a = j * var_c_res / var_w_res * env["delta"]
    theta_c = j * var_c_res / var_w * env["delta"]
    wc = j * cov_cz
    r2w = wc @ np.linalg.solve(variance_z, wc) / var_w
    source_slope = np.linalg.solve(variance_z, cov_cz)
    source_coef_c = np.r_[pi - mean_z @ source_slope, source_slope]
    h = np.concatenate(env["h"])
    mean_h = weights @ h
    cov_zh = (z - mean_z).T @ (weights * (h - mean_h))
    cov_wh = j * (weights * (c - pi)) @ h
    h_on_w = (cov_wh - wc @ np.linalg.solve(variance_z, cov_zh)) / var_w_res
    return {"theta_a": theta_a, "theta_c": theta_c, "r2_c": 1 - var_c_res / variance_c,
            "r2_w": r2w, "source_coef_c": source_coef_c, "h_on_w": h_on_w}


def generate(env, n, pi, size, rng, empirical=False):
    n1 = math.ceil(pi * n)
    c = np.r_[np.ones(n1), np.zeros(n - n1)]
    z = np.empty((size, n, 5))
    noise = np.empty((size, n, 6))
    for g, selection in ((1, slice(0, n1)), (0, slice(n1, n))):
        count = n1 if g else n - n1
        z[:, selection] = env["z"][g][rng.integers(len(env["z"][g]), size=(size, count))]
        if empirical:
            noise[:, selection] = env["empirical"][g][rng.integers(len(env["empirical"][g]), size=(size, count))]
        else:
            noise[:, selection] = rng.standard_normal((size, count, 6)) @ env["factor"][g].T
    k = np.concatenate((np.ones((size, n, 1)), z), axis=2)
    h = (z[:, :, 0] * z[:, :, 1] - k @ env["hcoef"]) / env["hsd"]
    y = c[None, :, None] * env["delta"] + z @ env["gamma"] + noise
    if env["interaction"]:
        y += .25 * h[:, :, None]
    return c, z, h, y


def projection(x, y):
    gram = np.swapaxes(x, 1, 2) @ x
    eigenvalues = np.linalg.eigvalsh(gram)
    valid = eigenvalues[:, 0] > 1e-10 * eigenvalues[:, -1]
    safe = gram.copy()
    safe[~valid] = np.eye(gram.shape[-1])
    inverse = np.linalg.inv(safe)
    weight = x @ inverse
    coefficients = np.swapaxes(weight, 1, 2) @ y
    residual = y - x @ coefficients
    leverage = np.sum(weight * x, axis=2)
    valid &= np.max(leverage, axis=1) < 1 - 1e-8
    return coefficients, residual, weight, leverage, valid


def hc3(x, y):
    coefficients, residual, weight, leverage, valid = projection(x, y)
    # HC3 uses the complete design leverage and N - rank(X) t reference.
    se = np.sqrt(np.sum((residual / np.maximum(1 - leverage[:, :, None], 1e-12) * weight[:, :, 1, None]) ** 2, axis=1))
    valid &= np.all(np.isfinite(se) & (se > 1e-12), axis=1)
    return coefficients[:, 1], se, valid, weight


def bh(p):
    order = np.argsort(p, axis=-1)
    ranked = np.take_along_axis(p, order, axis=-1)
    adjusted = ranked * p.shape[-1] / np.arange(1, p.shape[-1] + 1)
    adjusted = np.minimum.accumulate(adjusted[..., ::-1], axis=-1)[..., ::-1].clip(0, 1)
    return np.take_along_axis(adjusted, np.argsort(order, axis=-1), axis=-1)


def simulate(s):
    start_time = time.perf_counter()
    env = model_environment(s)
    n, repetitions = s["n"], s["repetitions"]
    pi = math.ceil(s["pi"] * n) / n
    pop = population(env, pi, s["se"], s["sp"], env["interaction"])
    pop_z = population(env, pi, s["se"], s["sp"], False)
    source_pop = population(env, .125, 1., 1., False)
    source_coefficient = np.vstack((np.zeros(6), env["gamma"])) + source_pop["source_coef_c"][:, None] * env["delta"]
    beta = np.empty((repetitions, 5, 6))
    se = np.empty_like(beta)
    divisor = np.empty_like(beta)
    target = np.empty_like(beta)
    valid = np.empty((repetitions, 5), dtype=bool)
    r2_observed = np.empty(repetitions)
    max_fwl_error = 0.
    for block, lo in enumerate(range(0, repetitions, BATCH)):
        hi = min(lo + BATCH, repetitions)
        b = hi - lo
        rng = stream(key(s), block)
        c, z, h, y = generate(env, n, s["pi"], b, rng, s["extension"] == "empirical_errors")
        probabilities = np.where(c == 1, s["se"], 1 - s["sp"])
        w = (rng.random((b, n)) < probabilities).astype(float)
        cov = np.concatenate((z, h[:, :, None]), axis=2) if env["interaction"] else z
        k = np.concatenate((np.ones((b, n, 1)), cov), axis=2)
        x = np.concatenate((np.ones((b, n, 1)), w[:, :, None], cov), axis=2)
        xo = x.copy()
        xo[:, :, 1] = c
        bo, so, vo, _ = hc3(xo, y)
        ba, sa, va, weights = hc3(x, y)
        _, yr, _, _, vk = projection(k, y)
        bb = np.einsum("bn,bnj->bj", weights[:, :, 1], yr)
        if np.any(va & vk):
            error = np.max(np.abs(bb[va & vk] - ba[va & vk]))
            max_fwl_error = max(max_fwl_error, float(error))
            if error > 1e-9:
                raise ArithmeticError("FWL coefficient identity failed")
        bc, sc, vc, _ = hc3(x[:, :, :2], yr)
        _, wr, _, _, vw = projection(k, w[:, :, None])
        denominator = np.sum((w - w.mean(1, keepdims=True)) ** 2, axis=1)
        r2_observed[lo:hi] = 1 - np.sum(wr[:, :, 0] ** 2, axis=1) / np.maximum(denominator, 1e-30)
        if env["interaction"]:
            _, zs, hs, ys = generate(env, 1000, .125, b, rng)
            ks = np.concatenate((np.ones((b, 1000, 1)), zs, hs[:, :, None]), axis=2)
            source_fit, _, _, _, vs = projection(ks, ys)
            yd = y - k @ source_fit
            xd = np.concatenate((np.ones((b, n, 1)), w[:, :, None], z), axis=2)
            bd, sd, vd, _ = hc3(xd, yd)
            vd &= vs
            td = pop_z["theta_a"][None, :] + (.25 - source_fit[:, -1]) * pop_z["h_on_w"]
            restored = np.einsum("bn,bnj->bj", weights[:, :, 1], yd)
            if np.any(va & vs):
                err = np.max(np.abs(restored[va & vs] - ba[va & vs]))
                max_fwl_error = max(max_fwl_error, float(err))
                if err > 1e-9:
                    raise ArithmeticError("Matched source-space identity failed")
        else:
            yd = y - k @ source_coefficient
            bd = np.einsum("bn,bnj->bj", weights[:, :, 1], yd)
            sd, vd = sa, va
            td = np.broadcast_to(pop["theta_a"], (b, 6))
        beta[lo:hi] = np.stack((bo, ba, bb, bc, bd), axis=1)
        se[lo:hi] = np.stack((so, sa, sa, sc, sd), axis=1)
        divisor[lo:hi] = np.stack((y.std(1, ddof=1), y.std(1, ddof=1), yr.std(1, ddof=1), yr.std(1, ddof=1), yd.std(1, ddof=1)), axis=1)
        target[lo:hi] = np.stack((np.broadcast_to(env["delta"], (b, 6)), np.broadcast_to(pop["theta_a"], (b, 6)),
                                np.broadcast_to(pop["theta_a"], (b, 6)), np.broadcast_to(pop["theta_c"], (b, 6)), td), axis=1)
        valid[lo:hi] = np.stack((vo, va, va & vk, vc & vk & vw, vd), axis=1)
    full_rank = 8 if env["interaction"] else 7
    dfs = np.array([n - full_rank, n - full_rank, n - full_rank, n - 2, n - 7])
    pvalues = 2 * stats.t.sf(np.abs(beta / se), dfs[None, :, None])
    pvalues[~valid] = 1.
    return {"beta": beta, "se": se, "divisor": divisor, "projection_target": target, "valid": valid,
            "p": pvalues, "df": dfs, "r2_observed": r2_observed,
            "max_fwl_error": np.array(max_fwl_error), "elapsed_seconds": np.array(time.perf_counter() - start_time)}


def expectation(values):
    values = np.asarray(values, dtype=float)
    mean = np.mean(values)
    mcse = np.std(values, ddof=1) / math.sqrt(len(values)) if len(values) > 1 else float("nan")
    return float(mean), float(mcse)


def summaries(s, data):
    rows, panels = [], []
    delta = np.asarray(EFFECTS[s["effect"]])
    q = bh(data["p"])
    base = {**s, "scenario": key(s), "design": key(s, True)}
    for m, method in enumerate(METHODS):
        ok = data["valid"][:, m]
        if not np.any(ok):
            raise ArithmeticError(f"No estimable repetitions: {key(s)} / {method}")
        rejected = q[:, m] < .05
        number = rejected.sum(1)
        false = rejected[:, delta == 0].sum(1)
        fdp = false / np.maximum(number, 1)
        panel = {**base, "method": method, "n_valid": int(ok.sum()), "n_failed": int((~ok).sum())}
        for name, values in {"fdr": fdp, "fwer": false > 0, "any_discovery": number > 0,
                             "true_discovery_fraction": rejected[:, delta != 0].mean(1) if np.any(delta != 0) else np.zeros(len(ok))}.items():
            panel[name], panel[name + "_mcse"] = expectation(values)
        panels.append(panel)
        for scale in ("native", "subgroup_standardised"):
            div = np.ones_like(data["divisor"][ok, m]) if scale == "native" else data["divisor"][ok, m]
            b, se = data["beta"][ok, m] / div, data["se"][ok, m] / div
            truth = delta / div
            projection_truth = data["projection_target"][ok, m] / div
            t95, t90 = stats.t.ppf([.975, .95], data["df"][m])
            for j, protein in enumerate(PROTEINS):
                row = {**base, "method": method, "scale": scale, "protein": protein,
                       "true_effect_reference": float(delta[j]), "n_valid": int(ok.sum()),
                       "n_failed": int((~ok).sum()), "failure_rate": float(np.mean(~ok))}
                error, projection_error = b[:, j] - truth[:, j], b[:, j] - projection_truth[:, j]
                metrics = {"mean_beta": b[:, j], "mean_scientific_target": truth[:, j],
                           "mean_projection_target": projection_truth[:, j], "bias_target": error,
                           "bias_projection": projection_error, "mean_se": se[:, j],
                           "mean_divisor": div[:, j], "coverage_95_target": abs(error) <= t95 * se[:, j],
                           "coverage_95_projection": abs(projection_error) <= t95 * se[:, j],
                           "width_95": 2*t95*se[:, j], "discovery_nominal": data["p"][:, m, j] < .05,
                           "discovery_bh6": q[:, m, j] < .05}
                for margin in (.3, .5):
                    name = str(round(margin * 100)).zfill(3)
                    translated = margin / div[:, j]
                    for label, decision in (("exclusion_", b[:, j] + t95 * se[:, j] < translated),
                                            ("equivalence_", abs(b[:, j]) + t90 * se[:, j] < translated),
                                            ("untranslated_equivalence_", abs(b[:, j]) + t90 * se[:, j] < margin)):
                        decisions = np.zeros(len(ok), dtype=bool)
                        decisions[ok] = decision
                        metrics[label + name] = decisions
                for name, values in metrics.items():
                    row[name], row[name + "_mcse"] = expectation(values)
                row["rmse_target"] = float(np.sqrt(np.mean(error ** 2)))
                row["rmse_target_mcse"] = float(np.std(error ** 2, ddof=1) / math.sqrt(len(error)) / (2*row["rmse_target"])) if row["rmse_target"] else 0.
                row["empirical_se"] = float(np.std(b[:, j], ddof=1))
                centred = b[:, j] - b[:, j].mean()
                row["empirical_se_mcse"] = float(np.std(centred ** 2, ddof=1) / math.sqrt(len(centred)) / (2*row["empirical_se"]))
                row["se_ratio"] = row["mean_se"] / row["empirical_se"]
                row["relative_bias_target"] = row["bias_target"] / row["mean_scientific_target"] if delta[j] else None
                rows.append(row)
    return pd.DataFrame(rows), pd.DataFrame(panels)


def run_cell(s, directory, signature):
    directory = Path(directory)
    dest = directory / "cells" / key(s)
    dest.mkdir(parents=True, exist_ok=True)
    complete = dest / "complete.json"
    if complete.exists():
        record = json.loads(complete.read_text())
        if record["signature"] != signature or record["repetitions"] != s["repetitions"]:
            raise ValueError(f"Incompatible checkpoint: {dest}")
        return {"scenario": key(s), "elapsed_seconds": record["elapsed_seconds"], "cached": True}
    data = simulate(s)
    estimates, panels = summaries(s, data)
    estimates.to_csv(dest / "estimates.csv", index=False)
    panels.to_csv(dest / "panels.csv", index=False)
    np.savez_compressed(dest / "replicates.npz", **data)
    record = {**s, "signature": signature, "scenario": key(s),
              "elapsed_seconds": float(data["elapsed_seconds"]), "max_fwl_error": float(data["max_fwl_error"]),
              "failed_by_method": dict(zip(METHODS, map(int, (~data["valid"]).sum(0))))}
    write_json(complete, record)
    return record


@contextlib.contextmanager
def output_lock(directory):
    import fcntl
    Path(directory).mkdir(parents=True, exist_ok=True)
    with (Path(directory) / ".run.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield


def execute(cells, calibration, directory, workers, signature):
    started = time.perf_counter()
    with ProcessPoolExecutor(max_workers=workers, initializer=initialise, initargs=(str(calibration),)) as executor:
        jobs = {executor.submit(run_cell, s, str(directory), signature): s for s in cells}
        for count, future in enumerate(as_completed(jobs), 1):
            result = future.result()
            record = {"status": "RUNNING", "completed_in_batch": count, "batch_total": len(cells),
                      "elapsed_seconds": time.perf_counter() - started, "last_completed": result["scenario"]}
            write_json(Path(directory) / "status.json", record)
            print(json.dumps(record), flush=True)


def aggregate(directory, cells):
    directory = Path(directory)
    estimate_frames, panel_frames, paired_rows = [], [], []
    group = {}
    for s in cells:
        folder = directory / "cells" / key(s)
        if not (folder / "complete.json").exists():
            raise FileNotFoundError(f"Incomplete cell: {key(s)}")
        estimate_frames.append(pd.read_csv(folder / "estimates.csv", keep_default_na=False, na_values=[""]))
        panel_frames.append(pd.read_csv(folder / "panels.csv", keep_default_na=False, na_values=[""]))
        group.setdefault(key(s, True), {})[s["cohort"]] = s
    estimates = pd.concat(estimate_frames, ignore_index=True)
    panels = pd.concat(panel_frames, ignore_index=True)
    estimates.to_csv(directory / "simulation_estimates.csv", index=False)
    panels.to_csv(directory / "simulation_panels.csv", index=False)
    for design, pair in group.items():
        if set(pair) != set(COHORTS):
            continue
        arrays = [np.load(directory / "cells" / key(pair[c]) / "replicates.npz") for c in COHORTS]
        pv = np.concatenate([a["p"] for a in arrays], axis=2)
        rejected = bh(pv) < .05
        truth = np.tile(np.asarray(EFFECTS[pair["INCOV"]["effect"]]), 2)
        for m, method in enumerate(METHODS):
            discoveries = rejected[:, m]
            false = discoveries[:, truth == 0].sum(1)
            fdp = false / np.maximum(discoveries.sum(1), 1)
            row = {"design": design, "method": method, "n_per_cohort": pair["INCOV"]["n"],
                   "se": pair["INCOV"]["se"], "sp": pair["INCOV"]["sp"], "effect": pair["INCOV"]["effect"],
                   "extension": pair["INCOV"]["extension"], "family_size": 12}
            for name, v in {"fdr": fdp, "fwer": false > 0,
                            "true_discovery_fraction": discoveries[:, truth != 0].mean(1) if np.any(truth != 0) else np.zeros(len(fdp))}.items():
                row[name], row[name + "_mcse"] = expectation(v)
            for j, protein in enumerate(PROTEINS):
                for ci, cohort in enumerate(COHORTS):
                    name = f"{cohort}_{protein}_power"
                    row[name], row[name + "_mcse"] = expectation(discoveries[:, ci*6+j])
            paired_rows.append(row)
        for a in arrays:
            a.close()
    pd.DataFrame(paired_rows).to_csv(directory / "simulation_bh12.csv", index=False)
    failures = panels.groupby("method").n_failed.sum().to_dict()
    record = {"status": "COMPLETE", "cells": len(cells), "dataset_repetitions": sum(s["repetitions"] for s in cells),
              "seed": SEED, "batch_size": BATCH, "failed_by_method": {k: int(v) for k, v in failures.items()},
              "python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform()}
    write_json(directory / "summary.json", record)
    write_json(directory / "status.json", record)
    return estimates


def tost_power(n, margin=.3):
    n1 = math.ceil(n * .25)
    ss = math.sqrt(1/n1 + 1/(n-n1))
    df = n - 2
    critical = stats.t.ppf(.95, df)
    upper = df * (margin / (ss * critical)) ** 2
    def integrand(v):
        return max(2*stats.norm.cdf(margin/ss-critical*math.sqrt(v/df))-1, 0) * stats.chi2.pdf(v, df)
    return integrate.quad(integrand, 0, upper, epsabs=1e-10)[0]


def validate(calibration, directory):
    initialise(calibration)
    checks = {}
    assert len(grid()) == 312
    test = scenario("INCOV", 100, .8, .85, "dense05", repetitions=128)
    env = model_environment(test)
    c, z, h, y = generate(env, 100, .25, 3, stream("validation_matrix"))
    w = np.broadcast_to(c, (3, 100)).copy()
    x = np.concatenate((np.ones((3, 100, 1)), w[:, :, None], z), 2)
    b, se, valid, _ = hc3(x, y)
    svd_errors = []
    for i in range(3):
        pinv = np.linalg.pinv(x[i])
        coefficients = pinv @ y[i]
        e = y[i] - x[i] @ coefficients
        leverage = np.sum(x[i] * pinv.T, axis=1)
        for j in range(6):
            sandwich = (pinv * (e[:, j]/(1-leverage)) ** 2) @ pinv.T
            svd_errors.append(abs(math.sqrt(sandwich[1, 1])-se[i, j]))
    assert max(svd_errors) < 1e-10
    checks["hc3_independent_svd_max_error"] = max(svd_errors)
    first = simulate(test)
    repeat = simulate(test)
    assert np.array_equal(first["beta"], repeat["beta"])
    with ProcessPoolExecutor(max_workers=2, initializer=initialise, initargs=(str(calibration),)) as pool:
        parallel = pool.submit(simulate, test).result()
    assert np.array_equal(first["beta"], parallel["beta"])
    checks["serial_parallel_identical"] = True
    checks["fwl_max_error"] = float(first["max_fwl_error"])
    ext = simulate(scenario("MYLC", 100, .8, .95, "dense03", "source_interaction", 128))
    checks["source_matched_space_max_error"] = float(ext["max_fwl_error"])
    assert checks["source_matched_space_max_error"] < 1e-9
    estimates, _ = summaries(test, first)
    for metric in ("discovery_nominal", "coverage_95_target", "equivalence_030", "exclusion_030"):
        assert np.allclose(estimates[estimates.scale == "native"][metric],
                           estimates[estimates.scale == "subgroup_standardised"][metric])
    checks["translated_scale_decisions_identical"] = True
    rng = stream("validation_normal_reference")
    for n, delta in ((100, 0.), (466, .3), (509, 0.)):
        n1, repetitions = math.ceil(.25*n), 10000
        difference = rng.normal(delta, math.sqrt(1/n1+1/(n-n1)), repetitions)
        variance = rng.chisquare(n-2, repetitions)/(n-2)
        se0 = np.sqrt(variance*(1/n1+1/(n-n1)))
        if n == 509:
            achieved = float(np.mean(abs(difference)+stats.t.ppf(.95,n-2)*se0 < .3))
            expected = tost_power(n)
        else:
            achieved = float(np.mean(abs(difference/se0)>stats.t.ppf(.975,n-2)))
            ncp = delta/math.sqrt(1/n1+1/(n-n1))
            threshold = stats.t.ppf(.975,n-2)
            expected = float(stats.nct.cdf(-threshold,n-2,ncp)+stats.nct.sf(threshold,n-2,ncp))
        tolerance = 5*math.sqrt(expected*(1-expected)/repetitions)
        assert abs(achieved-expected)<tolerance
        checks[f"normal_reference_n{n}"] = {"achieved":achieved,"analytic":expected,"mcse":tolerance/5}
    null = simulate(scenario("INCOV", 100, 1., 1., "null", repetitions=2000))
    null_summary, null_panel = summaries(scenario("INCOV", 100, 1., 1., "null", repetitions=2000), null)
    null_summary.to_csv(Path(directory)/"validation_null_estimates.csv", index=False)
    checks["calibrated_null_failures"] = int((~null["valid"]).sum())
    assert checks["calibrated_null_failures"] == 0
    checks["status"] = "PASS"
    checks["source_sha256"] = sha(__file__)
    checks["calibration_sha256"] = sha(calibration)
    write_json(Path(directory)/"validation.json", checks)
    print(json.dumps(checks, indent=2))



PROFILES = ((1., 1.), (.8, .95), (.8, .85))
MAX_N = 5000
BUDGET = 96


def strata():
    return list(itertools.product(COHORTS, PROFILES, ("null", "dense03", "dense05")))


def metrics(effect):
    if effect == "null":
        return ("exclusion_030", "exclusion_050", "equivalence_030", "equivalence_050")
    return ("discovery_nominal", "discovery_bh6")


def select_rows(frame, stratum):
    cohort, (se, sp), effect = stratum
    return frame[(frame.cohort == cohort) & (frame.se == se) & (frame.sp == sp) & (frame.effect == effect)]


def probability_curve(frame, protein, metric):
    d = frame[frame.protein == protein].sort_values(["n", "repetitions"]).drop_duplicates("n", keep="last")
    x = d.n.to_numpy(float)
    p = d[metric].to_numpy(float)
    e = d[metric + "_mcse"].to_numpy(float)
    return x, p, e


def bracket(x, p, e):
    lower = x[p + 1.96 * e < .8]
    upper = x[p - 1.96 * e >= .8]
    lo = float(lower.max()) if len(lower) else 0.
    hi = float(upper.min()) if len(upper) else float("inf")
    return lo, hi


def proposal(frame, stratum):
    d = select_rows(frame, stratum)
    existing = set(d.n)
    options = []
    for protein, metric in itertools.product(PROTEINS, metrics(stratum[2])):
        x, p, e = probability_curve(d, protein, metric)
        lo, hi = bracket(x, p, e)
        if hi <= lo:
            continue
        if math.isinf(hi):
            if x.max() >= MAX_N:
                continue
            n = min(MAX_N, int(math.ceil(x.max() * 1.6 / 10) * 10))
            score = 2.0
        elif lo == 0:
            if hi <= 40:
                continue
            n = max(40, int(round(hi / 1.6 / 10) * 10))
            score = 1.5
        else:
            if hi - lo <= 10:
                continue
            n = int(round(math.sqrt(lo * hi) / 10) * 10)
            n = min(int(hi - 1), max(int(lo + 1), n))
            score = math.log(hi / lo)
        if n in existing:
            alternatives = [v for v in range(max(40, int(lo) + 10), min(MAX_N, int(hi) if math.isfinite(hi) else MAX_N) + 1, 10) if v not in existing]
            if not alternatives:
                continue
            n = min(alternatives, key=lambda v: abs(v - n))
        options.append((score, n))
    return max(options) if options else None


def load_estimates(directory, primary):
    frames = [primary]
    for path in sorted((directory / "cells").glob("*/estimates.csv")):
        d = pd.read_csv(path, keep_default_na=False, na_values=[""])
        frames.append(d[(d.method == "a_joint") & (d.scale == "native")])
    return pd.concat(frames, ignore_index=True)


def recruits(n):
    if not math.isfinite(n):
        return None
    n = int(n)
    positive = math.ceil(.25 * n)
    return math.ceil(positive / .8) + math.ceil((n - positive) / .8)


def interpolate(x, p):
    # Monotone interpolation supplies a planning estimate; Monte Carlo brackets are retained.
    p = np.maximum.accumulate(p)
    if p[0] >= .8:
        return float(x[0]), "at_or_below_smallest_evaluated_N"
    if p[-1] < .8:
        return None, "above_largest_evaluated_N"
    k = int(np.flatnonzero(p >= .8)[0])
    value = x[k-1] + (.8 - p[k-1]) * (x[k] - x[k-1]) / (p[k] - p[k-1])
    return float(math.ceil(value / 10) * 10), "interpolated_10_participant_grid"


def tabulate(frame, directory):
    rows = []
    for stratum in strata():
        cohort, (se, sp), effect = stratum
        d = select_rows(frame, stratum)
        for protein, metric in itertools.product(PROTEINS, metrics(effect)):
            x, p, e = probability_curve(d, protein, metric)
            lo, hi = bracket(x, p, e)
            n, status = interpolate(x, p)
            baseline = 466 if metric.endswith("030") or effect == "dense03" else 169
            if metric.startswith("equivalence"):
                baseline = 509 if metric.endswith("030") else 185
            if metric == "discovery_bh6":
                comparison = "uncorrected_two_group_detection_reference"
            else:
                comparison = "two_group_analytic_reference"
            if effect == "null" and (se, sp) != (1., 1.):
                scope = "proxy_endpoint_only; latent_effect_error_control_not_established"
            elif effect != "null" and (se, sp) != (1., 1.):
                scope = "discovery_of_attenuated_proxy_association; not_latent_effect_coverage"
            else:
                scope = "correctly_measured_phenotype_under_calibrated_model"
            rows.append({"cohort": cohort, "se": se, "sp": sp, "effect": effect,
                         "protein": protein, "objective": metric, "target_probability": .80,
                         "estimated_N80": n, "estimate_type": status,
                         "N_below_80_MC95": lo if lo else None,
                         "N_above_80_MC95": hi if math.isfinite(hi) else None,
                         "maximum_evaluated_N": int(x.max()), "recruits_20pct_attrition": recruits(n) if n else None,
                         "analytic_reference_N": baseline, "reference_type": comparison,
                         "conditional_N_ratio": n / baseline if n else None, "interpretation": scope})
    table = pd.DataFrame(rows)
    table.to_csv(directory / "Table_2_conditional_design.csv", index=False)
    frame.to_csv(directory / "conditional_design_points.csv", index=False)
    write_json(directory / "conditional_design_summary.json", {
        "additional_cells": len(list((directory/"cells").glob("*/complete.json"))),
        "repetitions_per_additional_cell": 20000, "maximum_N": MAX_N,
        "interval_definition": "N with upper Monte Carlo 95% limit below 0.8; N with lower limit at least 0.8",
        "rows": json.loads(table.to_json(orient="records"))})


def design_search(workers):
    calibration = ROOT / "data/ademp_calibration.json"
    directory = ROOT / "output/design"
    directory.mkdir(parents=True, exist_ok=True)
    primary = pd.read_csv(ROOT / "output/ademp/simulation_estimates.csv", keep_default_na=False, na_values=[""])
    primary = primary[(primary.extension == "main") & (primary.method == "a_joint") & (primary.scale == "native")]
    signature = hashlib.sha256((sha(__file__) + sha(calibration)).encode()).hexdigest()
    with output_lock(directory):
        while True:
            frame = load_estimates(directory, primary)
            completed = len(list((directory / "cells").glob("*/complete.json")))
            remaining = BUDGET - completed
            if remaining <= 0:
                break
            candidates = []
            for stratum in strata():
                item = proposal(frame, stratum)
                if item:
                    priority, n = item
                    cohort, (se, sp), effect = stratum
                    candidates.append((priority, scenario(cohort, n, se, sp, effect, "design", 20000)))
            if not candidates:
                break
            candidates.sort(key=lambda item: (-item[0], key(item[1])))
            cells = [s for _, s in candidates[:remaining]]
            execute(cells, calibration, directory, workers, signature)
        frame = load_estimates(directory, primary)
        tabulate(frame, directory)
        write_json(directory / "status.json", {"status": "COMPLETE", "additional_cells": len(list((directory/"cells").glob("*/complete.json"))), "budget": BUDGET})



def summarise_all():
    folder = ROOT / "output/ademp"
    plan = json.loads((folder / "design.json").read_text())
    paired, directional, key_findings = [], [], []
    for cell in plan["cells"]:
        with np.load(folder / "cells" / key(cell) / "replicates.npz") as d:
            truth = np.asarray(EFFECTS[cell["effect"]])
            b, se, valid = d["beta"], d["se"], d["valid"]
            discoveries = bh(d["p"]) < .05
            t95 = stats.t.ppf(.975, d["df"])[None, :, None]
            t90 = stats.t.ppf(.95, d["df"])[None, :, None]
            rates = {"coverage_target": abs(b-truth) <= t95*se,
                     "discovery_bh6": discoveries,
                     "exclusion_030": b+t95*se < .3,
                     "equivalence_030": abs(b)+t90*se < .3}
            if cell["extension"] == "main" and cell["effect"] in ("dense03", "dense05") and cell["se"] == .8 and cell["sp"] == .85 and cell["n"] in (100, 1000):
                outcomes = {**rates, "coverage_projection": abs(b-d["projection_target"]) <= t95*se,
                            "exclusion_050": b+t95*se < .5, "equivalence_050": abs(b)+t90*se < .5}
                for metric, values in outcomes.items():
                    array=values[:,1].astype(float)
                    mean,mcse=expectation(array.mean(1))
                    key_findings.append({**cell,"metric":metric,"six_protein_mean":mean,"mcse":mcse,
                                         "protein_minimum":float(array.mean(0).min()),"protein_maximum":float(array.mean(0).max())})
            for m, method in enumerate(METHODS):
                for j, protein in enumerate(PROTEINS):
                    correct = discoveries[:, m, j] & (b[:, m, j] * truth[j] > 0)
                    p, e = expectation(correct)
                    directional.append({**cell, "method": method, "protein": protein,
                                        "directional_discovery_bh6": p, "mcse": e})
                    if method == "a_joint":
                        continue
                    ok = valid[:, m] & valid[:, 1]
                    for metric, values in rates.items():
                        diff = values[ok, m, j].astype(float) - values[ok, 1, j].astype(float)
                        mean, mcse = expectation(diff)
                        paired.append({**cell, "method": method, "reference_method": "a_joint",
                                       "protein": protein, "metric": metric, "paired_difference": mean,
                                       "mcse": mcse, "paired_valid_repetitions": int(ok.sum())})
    pd.DataFrame(paired).to_csv(folder / "paired_method_comparisons.csv", index=False)
    pd.DataFrame(directional).to_csv(folder / "directional_discovery.csv", index=False)
    pd.DataFrame(key_findings).to_csv(folder / "key_findings.csv", index=False)
    initialise(ROOT / "data/ademp_calibration.json")
    calibration_rows, error_rows = [], []
    for cohort in COHORTS:
        s = scenario(cohort, 1000, 1., 1., "dense05")
        env = model_environment(s)
        item = CALIBRATION["cohorts"][cohort]
        pooled = .75*env["covariance"][0] + .25*env["covariance"][1]
        correlation = pooled/np.sqrt(np.diag(pooled)[:, None]*np.diag(pooled)[None, :])
        pd.DataFrame(correlation, index=PROTEINS, columns=PROTEINS).to_csv(folder / f"{cohort}_residual_correlation.csv")
        for j, protein in enumerate(PROTEINS):
            calibration_rows.append({"cohort": cohort, "protein": protein,
                                     "residual_variance_ratio": env["covariance"][1][j,j]/env["covariance"][0][j,j],
                                     "nuisance_partial_R2": item["nuisance_partial_r2"][j],
                                     "cough_R2_Z": population(env,.25,1.,1.)["r2_c"]})
        for se, sp in ERROR_PROFILES:
            pop = population(env,.25,se,sp)
            pw = .25*se+.75*(1-sp)
            error_rows.append({"cohort": cohort, "sensitivity": se, "specificity": sp,
                               "latent_prevalence": .25, "observed_prevalence": pw,
                               "positive_predictive_value": .25*se/pw,
                               "joint_effect_ratio": pop["theta_a"][0]/.5,
                               "outcome_only_residual_effect_ratio": pop["theta_c"][0]/.5,
                               "R2_C_Z": pop["r2_c"], "R2_W_Z": pop["r2_w"]})
    pd.DataFrame(calibration_rows).to_csv(folder / "empirical_calibration.csv", index=False)
    pd.DataFrame(error_rows).to_csv(folder / "misclassification_projection.csv", index=False)
    design_folder=ROOT/"output/design"
    if (design_folder/"status.json").exists() and json.loads((design_folder/"status.json").read_text()).get("status")=="COMPLETE":
        table=pd.read_csv(design_folder/"Table_2_conditional_design.csv",keep_default_na=False,na_values=[""])
        columns=["residual_variance_ratio","nuisance_partial_R2","cough_R2_Z","N_ratio_to_perfect_phenotype"]
        table=table.drop(columns=[c for c in columns if c in table])
        table=table.merge(pd.DataFrame(calibration_rows),on=["cohort","protein"],validate="many_to_one")
        perfect=table[(table.se==1)&(table.sp==1)][["cohort","protein","effect","objective","estimated_N80"]].rename(columns={"estimated_N80":"perfect_N80"})
        table=table.merge(perfect,on=["cohort","protein","effect","objective"],validate="many_to_one")
        table["N_ratio_to_perfect_phenotype"]=table.estimated_N80/table.perfect_N80
        table=table.drop(columns="perfect_N80")
        table.to_csv(design_folder/"Table_2_conditional_design.csv",index=False)
        record=json.loads((design_folder/"conditional_design_summary.json").read_text())
        record["rows"]=json.loads(table.to_json(orient="records",double_precision=15))
        write_json(design_folder/"conditional_design_summary.json",record)
        ranges=table.groupby(["cohort","se","sp","effect","objective"],dropna=False).agg(
            N80_min=("estimated_N80","min"),N80_max=("estimated_N80","max"),
            phenotype_inflation_min=("N_ratio_to_perfect_phenotype","min"),
            phenotype_inflation_max=("N_ratio_to_perfect_phenotype","max"),
            resolved_proteins=("estimated_N80","count")).reset_index()
        ranges.to_csv(design_folder/"Table_2_design_ranges.csv",index=False)
    tables = {name: json.loads(pd.read_csv(folder/f"{name}.csv", keep_default_na=False, na_values=[""]).to_json(orient="records", double_precision=15))
              for name in ("simulation_estimates", "simulation_panels", "simulation_bh12",
                           "paired_method_comparisons", "directional_discovery",
                           "empirical_calibration", "misclassification_projection", "key_findings")}
    tables["run"] = json.loads((folder / "summary.json").read_text())
    if (design_folder/"status.json").exists() and json.loads((design_folder/"status.json").read_text()).get("status")=="COMPLETE":
        for name in ("estimates","panels"):
            frames=[pd.read_csv(p,keep_default_na=False,na_values=[""]) for p in sorted((design_folder/"cells").glob(f"*/{name}.csv"))]
            combined=pd.concat(frames,ignore_index=True)
            combined.to_csv(design_folder/f"conditional_simulation_{name}.csv",index=False)
            tables[f"conditional_simulation_{name}"]=json.loads(combined.to_json(orient="records",double_precision=15))
        tables["conditional_design"]=json.loads((design_folder/"conditional_design_summary.json").read_text())
    write_json(folder / "simulation_summary.json", tables)
    complete = [json.loads(p.read_text()) for p in (folder / "cells").glob("*/complete.json")]
    write_json(folder / "computational_checks.json", {
        "complete_cells": len(complete), "unique_cells": len({x["scenario"] for x in complete}),
        "dataset_repetitions": sum(x["repetitions"] for x in complete),
        "maximum_FWL_difference": max(x["max_fwl_error"] for x in complete),
        "estimated_core_seconds": sum(x["elapsed_seconds"] for x in complete),
        "failure_counts": {m:sum(x["failed_by_method"][m] for x in complete) for m in METHODS}})

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("command", choices=("prepare", "validate", "run", "design", "summarise"))
    parser.add_argument("--source", type=Path)
    parser.add_argument("--calibration", type=Path, default=ROOT/"data/ademp_calibration.json")
    parser.add_argument("--output", type=Path, default=ROOT/"output/ademp")
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--repetitions", type=int, default=10000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.command == "prepare":
        if args.source is None:
            parser.error("prepare requires --source")
        prepare(args.source, args.calibration)
    elif args.command == "validate":
        validate(args.calibration, args.output)
    elif args.command == "design":
        design_search(args.workers)
    elif args.command == "summarise":
        aggregate(args.output, grid(args.repetitions))
        summarise_all()
    else:
        cells = grid(args.repetitions)
        signature = hashlib.sha256((sha(__file__)+sha(args.calibration)+json.dumps(cells,sort_keys=True)).encode()).hexdigest()
        validation = json.loads((args.output/"validation.json").read_text())
        if validation["status"] != "PASS" or validation["source_sha256"] != sha(__file__) or validation["calibration_sha256"] != sha(args.calibration):
            raise ValueError("Validation must match the analysis code and calibration")
        with output_lock(args.output):
            write_json(args.output/"design.json", {"cells":cells,"signature":signature,"code_sha256":sha(__file__),"calibration_sha256":sha(args.calibration)})
            execute(cells,args.calibration,args.output,args.workers,signature)
            aggregate(args.output,cells)


if __name__ == "__main__":
    main()
