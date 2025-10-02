"""
Fix short, discordant local-ancestry phase errors between two phased haplotypes
and emit cleaned `.bed` files.

Overview
--------
Given two BED/TSV files with local ancestry tracts for haplotype 0 and
haplotype 1 of the **same individual**, this tool:

1) Standardizes inputs (computes `window_size = egpos - sgpos` if missing,
   sorts windows, and assigns per-chromosome `id` indices).
2) For each chromosome present in both haplotypes, finds *candidate phase
   switch points* where:
      • the windows have exactly the same coordinates (`chm`, `sgpos`, `egpos`),
      • the ancestry labels disagree between haps (`ancestry_0 != ancestry_1`),
      • and the window is short (default `threshold_cm = 2.0`).
3) For each candidate, evaluates whether **swapping** the ancestry labels
   between haplotypes would increase short-range continuity, using the score:
      cont = (l0 == mid0) + (r0 == mid0) + (l1 == mid1) + (r1 == mid1)
   The swap is applied **only if** the continuity score would strictly improve
   (`cont_switch > cont_noswitch`). Boundary windows (missing neighbors) are
   handled safely.
4) After switching, collapses consecutive windows with the same `ancestry`
   into single continuous segments (per chromosome).
5) Writes fixed haplotype files preserving the input column order but **dropping**
   `id` and `window_size` in the outputs.

Input format
------------
Two tab-separated files (hap0, hap1). Required columns:
    - chm      : chromosome label (e.g., "chr1", "chrX")
    - sgpos    : genetic start position (cM)
    - egpos    : genetic end position (cM)
    - ancestry : local ancestry label (e.g., "AFR", "EUR", "NAT")

Optional columns (kept and propagated when present):
    - spos, epos : physical start/end (bp)
    - window_size: (cM); computed if missing
    - id         : per-chromosome index; re-assigned internally

Assumptions & caveats
---------------------
- Inputs represent *the same* individual’s two phased haplotypes on the same
  coordinate system; candidate detection requires exact (`chm`, `sgpos`, `egpos`)
  matches between haplotypes.
- Only chromosomes present in **both** inputs are processed.
- Continuity scoring looks only at immediate left/right neighbors within each hap.
- Ties (no continuity improvement) are **not** switched.
- Default short-window threshold is 2.0 cM (changeable inside `fix_phase_errors`).

Outputs
-------
Two fixed `.bed` (TSV) files:
- Column order matches the input, except `id` and `window_size` are dropped.
- Consecutive same-ancestry windows are merged into single segments.

CLI usage
---------
    python fix_phase.py <hap0_input.bed> <hap1_input.bed> \
                        <hap0_output.fixed.bed> <hap1_output.fixed.bed>

Example
-------
    python fix_phase.py NA19700_0.bed NA19700_1.bed \
                        NA19700_0.fixed.bed NA19700_1.fixed.bed

Implementation notes
--------------------
Core functions (unchanged):
- find_candidate_switch_points: build candidate table (exact-match, discordant, short)
- neighbor_ancestry           : get left/right ancestry by per-chr id
- evaluate_switch             : compute continuity scores and switch decision
- make_switches               : apply beneficial switches and log changes
- merge_same_ancestry_runs_single_chr: collapse identical consecutive labels
- fix_phase_errors            : orchestrate per-chromosome pipeline
- save_fixed_bed_preserve_order: write outputs preserving column order
"""

import numpy as np
import pandas as pd
import argparse

# --- Exact alignment: same chm, sgpos, egpos ---


def find_candidate_switch_points(
    tracts_0: pd.DataFrame, tracts_1: pd.DataFrame, threshold_cm: float = 2.0
) -> pd.DataFrame:
    """
    Build a table of candidate switch points:
      - exact coordinate match: same [chm, sgpos, egpos]
      - discordant ancestry: ancestry_0 != ancestry_1
      - short window on hap0: window_size_0 <= threshold_cm

    Assumes both inputs have columns:
      chm, sgpos, egpos, ancestry  (and optionally id, window_size)

    Returns a DataFrame sorted by [chm, sgpos, egpos] with suffixed columns
    (e.g., ancestry_0/1, id_0/1, window_size_0/1).
    """
    required = {"chm", "sgpos", "egpos", "ancestry"}
    missing0 = required - set(tracts_0.columns)
    missing1 = required - set(tracts_1.columns)
    if missing0:
        raise KeyError(f"tracts_0 missing columns: {missing0}")
    if missing1:
        raise KeyError(f"tracts_1 missing columns: {missing1}")

    t0 = tracts_0.copy()
    t1 = tracts_1.copy()

    # ensure numeric coords and window_size
    for df in (t0, t1):
        df["sgpos"] = pd.to_numeric(df["sgpos"], errors="raise")
        df["egpos"] = pd.to_numeric(df["egpos"], errors="raise")
        if "window_size" not in df.columns:
            df["window_size"] = df["egpos"] - df["sgpos"]

    exact = t0.merge(
        t1, on=["chm", "sgpos", "egpos"], how="inner", suffixes=("_0", "_1")
    )

    # discordant ancestry & short windows
    exact = exact[exact["ancestry_0"] != exact["ancestry_1"]]
    candidates = exact[exact["window_size_0"] <= threshold_cm].copy()

    return candidates.sort_values(["chm", "sgpos", "egpos"]).reset_index(drop=True)


# what is the ancestry to the left in the first haplotype?
# we have to look at id - 1
# --- Helpers to look up left/right neighbors on a haplotype ---
def neighbor_ancestry(df_hap, id_val):
    """
    Assume df_hap is already filtered to a single chromosome and has columns: id, anc.
    Returns (left_anc, right_anc) for the given integer id.
    """
    # left
    if id_val > 0:
        left = df_hap.loc[df_hap["id"] == id_val - 1, "ancestry"]
        left_anc = left.iloc[0] if len(left) else np.nan
    else:
        left_anc = np.nan

    # right
    right = df_hap.loc[df_hap["id"] == id_val + 1, "ancestry"]
    right_anc = right.iloc[0] if len(right) else np.nan

    return left_anc, right_anc


def evaluate_switch(row: pd.Series, tr0: pd.DataFrame, tr1: pd.DataFrame) -> dict:
    """
    Evaluate whether swapping the mid labels between haplotypes (id_0 <-> id_1)
    increases continuity, defined as:
      cont = (l0 == mid0) + (r0 == mid0) + (l1 == mid1) + (r1 == mid1)

    Returns a dict with decision and scores.
    """
    id0 = int(row["id_0"])
    id1 = int(row["id_1"])

    # current labels pulled from the tract tables (avoid stale merged values)
    m0 = tr0.loc[tr0["id"] == id0, "ancestry"].iloc[0]
    m1 = tr1.loc[tr1["id"] == id1, "ancestry"].iloc[0]

    # neighbors
    l0, r0 = neighbor_ancestry(tr0, id0)
    l1, r1 = neighbor_ancestry(tr1, id1)

    # nan-safe equality
    def eq(a, b):
        return (not pd.isna(a)) and (not pd.isna(b)) and (a == b)

    # continuity with/without swap
    cont_switch = int(eq(l0, m1)) + int(eq(r0, m1)) + int(eq(l1, m0)) + int(eq(r1, m0))
    cont_noswitch = (
        int(eq(l0, m0)) + int(eq(r0, m0)) + int(eq(l1, m1)) + int(eq(r1, m1))
    )

    return {
        "id_0": id0,
        "id_1": id1,
        "left0": l0,
        "cur0": m0,
        "right0": r0,
        "left1": l1,
        "cur1": m1,
        "right1": r1,
        "new0": m1,
        "new1": m0,  # proposed swapped labels
        "cont_noswitch": cont_noswitch,
        "cont_switch": cont_switch,
        "should_switch": (cont_switch > cont_noswitch),
    }


def make_switches(
    tracts_0: pd.DataFrame,
    tracts_1: pd.DataFrame,
    switch_candidates: pd.DataFrame,
    verbose: bool = True,
):
    """
    Evaluate and apply ancestry switches on exact-match, short, discordant windows.

    Parameters
    ----------
    tracts_0, tracts_1 : DataFrames
        Must include columns: ['id', 'ancestry'] (already single-chromosome).
    switch_candidates : DataFrame
        Output from find_candidate_switch_points; must include ['id_0','id_1'].
        (Optional for logs: ['chm','sgpos','egpos','ancestry_0','ancestry_1','window_size_0'])
    verbose : bool
        Print per-switch info and a summary.

    Returns
    -------
    new_tr0, new_tr1, changes_df
      new_tr0/new_tr1 are copies with applied switches.
      changes_df logs each applied change with before/after and scores.
    """
    required_cols = {"id", "ancestry"}
    if not required_cols.issubset(tracts_0.columns):
        miss = required_cols - set(tracts_0.columns)
        raise KeyError(f"tracts_0 missing columns: {miss}")
    if not required_cols.issubset(tracts_1.columns):
        miss = required_cols - set(tracts_1.columns)
        raise KeyError(f"tracts_1 missing columns: {miss}")
    if not {"id_0", "id_1"}.issubset(switch_candidates.columns):
        raise KeyError("switch_candidates must include 'id_0' and 'id_1'")

    # Work on copies
    tr0 = tracts_0.copy()
    tr1 = tracts_1.copy()

    # Deterministic order (if coords exist)
    order_cols = [
        c for c in ["chm", "sgpos", "egpos"] if c in switch_candidates.columns
    ]
    cands = (
        switch_candidates.sort_values(order_cols).reset_index(drop=True)
        if order_cols
        else switch_candidates.copy()
    )

    changes = []
    if verbose:
        print(f"Evaluating {len(cands)} candidate switch points...")

    for i, row in cands.iterrows():
        ev = evaluate_switch(row, tr0, tr1)  # uses current state

        if ev["should_switch"]:
            id0, id1 = ev["id_0"], ev["id_1"]
            old0, old1 = ev["cur0"], ev["cur1"]
            new0, new1 = ev["new0"], ev["new1"]

            # apply swap
            tr0.loc[tr0["id"] == id0, "ancestry"] = new0
            tr1.loc[tr1["id"] == id1, "ancestry"] = new1

            # log
            log_row = {
                "idx": i,
                "id_0": id0,
                "id_1": id1,
                "old0": old0,
                "new0": new0,
                "old1": old1,
                "new1": new1,
                "cont_noswitch": ev["cont_noswitch"],
                "cont_switch": ev["cont_switch"],
            }
            for c in [
                "chm",
                "sgpos",
                "egpos",
                "window_size_0",
                "ancestry_0",
                "ancestry_1",
            ]:
                if c in cands.columns:
                    log_row[c] = row[c]
            changes.append(log_row)

            if verbose:
                coord_str = (
                    f" @ {row['chm']}:{row['sgpos']}-{row['egpos']}"
                    if {"chm", "sgpos", "egpos"}.issubset(cands.columns)
                    else ""
                )
                print(
                    f"- switch {id0}/{id1}{coord_str}: "
                    f"{old0}->{new0} (hap0), {old1}->{new1} (hap1) | "
                    f"{ev['cont_noswitch']}→{ev['cont_switch']}"
                )

    changes_df = pd.DataFrame(changes)
    sort_cols = [
        c for c in ["chm", "sgpos", "egpos", "id_0", "id_1"] if c in changes_df.columns
    ]
    if not changes_df.empty and sort_cols:
        changes_df = changes_df.sort_values(sort_cols).reset_index(drop=True)

    if verbose:
        print(f"Applied {len(changes_df)} / {len(cands)} switches.")

    return tr0, tr1, changes_df


def merge_same_ancestry_runs_single_chr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse consecutive rows with the same `ancestry` into a single segment.
    Assumes `df` contains only one chromosome and has columns:
    ['chm','spos','epos','ancestry','sgpos','egpos','window_size','id'].

    Returns a DataFrame with the same columns, where:
      - 'spos'/'sgpos' are taken from the first row of each run,
      - 'epos'/'egpos' are taken from the last row of each run,
      - 'window_size' is recomputed,
      - 'id' is reindexed (0..n-1) for this chromosome.
    """
    # Ensure order so "consecutive" is well-defined
    df = df.sort_values(["sgpos", "egpos"]).reset_index(drop=True).copy()

    # New run starts whenever ancestry changes (single chromosome assumption)
    run_break = df["ancestry"].ne(df["ancestry"].shift()).fillna(True)
    df["run_id"] = run_break.cumsum()

    # Collapse runs
    out = (
        df.groupby("run_id", as_index=False)
        .agg(
            {
                "chm": "first",
                "spos": "first",
                "epos": "last",
                "ancestry": "first",
                "sgpos": "first",
                "egpos": "last",
            }
        )
        .drop(columns=["run_id"])
    )

    # Recompute size and reindex ids
    out["window_size"] = out["egpos"] - out["sgpos"]
    out["id"] = range(len(out))

    # Return in your expected column order
    return out[
        ["chm", "spos", "epos", "ancestry", "sgpos", "egpos", "window_size", "id"]
    ]


def fix_phase_errors(
    tracts_0: pd.DataFrame,
    tracts_1: pd.DataFrame,
    threshold_cm: float = 2.0,
    verbose: bool = True,
):
    # Window size
    for df in (tracts_0, tracts_1):
        df["sgpos"] = pd.to_numeric(df["sgpos"])
        df["egpos"] = pd.to_numeric(df["egpos"])
        df["window_size"] = df["egpos"] - df["sgpos"]

    # --- Sort and add per-haplotype IDs (restart per chromosome) ---
    def sort_and_add_id(df):
        df = df.sort_values(["chm", "sgpos", "egpos"]).reset_index(drop=True).copy()
        df["id"] = df.groupby("chm", sort=False).cumcount().astype(int)
        return df

    tracts_0 = sort_and_add_id(tracts_0)
    tracts_1 = sort_and_add_id(tracts_1)

    # Iterate per chromosome and apply switches
    common_chroms = sorted(set(tracts_0["chm"]).intersection(tracts_1["chm"]))

    fixed_tracts_0 = []
    fixed_tracts_1 = []
    all_changes = []

    for ch in common_chroms:
        tr0_chr = tracts_0[tracts_0["chm"] == ch].copy()
        tr1_chr = tracts_1[tracts_1["chm"] == ch].copy()

        # find candidates on this chromosome
        switch_candidates = find_candidate_switch_points(
            tr0_chr, tr1_chr, threshold_cm=threshold_cm
        )

        # apply switches
        new_tr0_chr, new_tr1_chr, changes_df = make_switches(
            tr0_chr, tr1_chr, switch_candidates, verbose=verbose
        )

        # merge same-ancestry runs
        new_tr0_chr = merge_same_ancestry_runs_single_chr(new_tr0_chr)
        new_tr1_chr = merge_same_ancestry_runs_single_chr(new_tr1_chr)

        fixed_tracts_0.append(new_tr0_chr)
        fixed_tracts_1.append(new_tr1_chr)

        if not changes_df.empty:
            if "chm" not in changes_df.columns:
                changes_df["chm"] = ch
            all_changes.append(changes_df)

    # Reassemble whole-genome dataframes
    fixed_tracts_0 = (
        pd.concat(fixed_tracts_0, ignore_index=True)
        if fixed_tracts_0
        else tracts_0.copy()
    )
    fixed_tracts_1 = (
        pd.concat(fixed_tracts_1, ignore_index=True)
        if fixed_tracts_1
        else tracts_1.copy()
    )

    # Optional: tidy sort
    fixed_tracts_0 = fixed_tracts_0.sort_values(
        ["chm", "sgpos", "egpos", "id"]
    ).reset_index(drop=True)
    fixed_tracts_1 = fixed_tracts_1.sort_values(
        ["chm", "sgpos", "egpos", "id"]
    ).reset_index(drop=True)

    # Combined change log
    change_log = (
        pd.concat(all_changes, ignore_index=True) if all_changes else pd.DataFrame()
    )

    if verbose:
        print(f"\nTotal switches applied: {len(change_log)}")

    return fixed_tracts_0, fixed_tracts_1, change_log


def save_fixed_bed_preserve_order(
    df: pd.DataFrame, out_path: str, drop_cols=("id", "window_size")
):
    cols = [c for c in df.columns if c not in drop_cols]  # preserves original order
    df[cols].to_csv(
        out_path, sep="\t", index=False
    )  # add header=False if your pipeline needs no header


def main():
    parser = argparse.ArgumentParser(
        description="Fix short, discordant ancestry tract phase errors between two haplotypes and write fixed .bed files."
    )
    parser.add_argument("hap0_input", help="Input BED/TSV for haplotype 0")
    parser.add_argument("hap1_input", help="Input BED/TSV for haplotype 1")
    parser.add_argument("hap0_output", help="Output .bed path for fixed haplotype 0")
    parser.add_argument("hap1_output", help="Output .bed path for fixed haplotype 1")
    args = parser.parse_args()

    # Read inputs
    tracts_0 = pd.read_csv(args.hap0_input, sep="\t")
    tracts_1 = pd.read_csv(args.hap1_input, sep="\t")

    # Run fix with your defaults
    fixed_tracts_0, fixed_tracts_1, change_log = fix_phase_errors(
        tracts_0, tracts_1, threshold_cm=2.0, verbose=True
    )

    # Save outputs (preserve column order, drop id/window_size)
    save_fixed_bed_preserve_order(fixed_tracts_0, args.hap0_output)
    save_fixed_bed_preserve_order(fixed_tracts_1, args.hap1_output)

    print(f"Done. Wrote:\n  - {args.hap0_output}\n  - {args.hap1_output}")
    print(f"Total switches applied: {len(change_log)}")


if __name__ == "__main__":
    main()
