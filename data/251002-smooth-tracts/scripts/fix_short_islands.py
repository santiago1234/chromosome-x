#!/usr/bin/env python3
import argparse
import pandas as pd


def smooth_islands_chr(df: pd.DataFrame,
                       short_cm: float = 1.0,
                       edge_support_cm: float = 2.0) -> pd.DataFrame:
    """
    Smooth short tracts on a SINGLE chromosome:
      - Left edge:  if first < short_cm and next >= edge_support_cm → flip to next ancestry
      - Middle:     if tract < short_cm and both neighbors exist & match → flip to that ancestry
      - Right edge: if last < short_cm and prev >= edge_support_cm → flip to prev ancestry
    Returns a copy with possibly updated 'ancestry'.
    """
    if df.empty:
        return df.copy()

    d = df.sort_values(['sgpos', 'egpos']).reset_index(drop=False).copy()  # keep original order via 'index'
    if 'window_size' not in d.columns:
        d['window_size'] = d['egpos'] - d['sgpos']

    anc = d['ancestry'].to_numpy()
    ws  = d['window_size'].to_numpy()
    n   = len(d)
    new = anc.copy()

    # case 1: left edge
    if n >= 2 and ws[0] < short_cm and ws[1] >= edge_support_cm:
        new[0] = anc[1]

    # case 2: middle islands
    if n >= 3:
        left  = anc[:-2]
        mid   = anc[1:-1]
        right = anc[2:]
        short = ws[1:-1] < short_cm
        mask  = short & (left == right) & (mid != left)
        new[1:-1][mask] = left[mask]

    # case 3: right edge
    if n >= 2 and ws[-1] < short_cm and ws[-2] >= edge_support_cm:
        new[-1] = anc[-2]

    d['ancestry'] = new
    d = d.set_index('index').sort_index()  # restore original ordering
    return d


def merge_runs_chr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge consecutive rows with the same 'ancestry' on a SINGLE chromosome.
    Keeps first 'spos/sgpos' and last 'epos/egpos', recomputes 'window_size',
    and reindexes 'id' from 0..n-1.
    """
    df = df.sort_values(['sgpos', 'egpos']).reset_index(drop=True).copy()

    run_break = df['ancestry'].ne(df['ancestry'].shift()).fillna(True)
    df['run_id'] = run_break.cumsum()

    out = (
        df.groupby('run_id', as_index=False)
          .agg({
              'chm': 'first',
              'spos': 'first',
              'epos': 'last',
              'ancestry': 'first',
              'sgpos': 'first',
              'egpos': 'last',
          })
          .drop(columns=['run_id'])
    )

    out['window_size'] = out['egpos'] - out['sgpos']
    out['id'] = range(len(out))
    return out[['chm','spos','epos','ancestry','sgpos','egpos','window_size','id']]


def smooth_and_merge_chr(tracts: pd.DataFrame,
                         short_cm: float = 1.0,
                         edge_support_cm: float = 2.0,
                         preview_n: int = 5) -> pd.DataFrame:
    """
    For a SINGLE chromosome: smooth short islands, print a brief summary,
    then merge consecutive same-ancestry runs. Returns the merged dataframe.
    """
    tr = tracts.copy()
    if 'window_size' not in tr.columns:
        tr['window_size'] = tr['egpos'] - tr['sgpos']

    # smooth
    tr_sm = smooth_islands_chr(tr, short_cm=short_cm, edge_support_cm=edge_support_cm)

    # summary
    changed = tr_sm['ancestry'].to_numpy() != tr['ancestry'].to_numpy()
    n_changed = int(changed.sum())
    print(f"short-tract flips applied: {n_changed} / {len(tr)} ({changed.mean():.1%})")
    if n_changed and preview_n > 0:
        prev = tr.loc[changed, ['sgpos','egpos','ancestry']].copy()
        prev['new_ancestry'] = tr_sm.loc[changed, 'ancestry'].values
        print(prev.head(preview_n))

    # merge
    merged = merge_runs_chr(tr_sm)
    print(f"before: {len(tr)} tracts | after smooth: {len(tr_sm)} | after merge: {len(merged)}")
    return merged


def smooth_and_merge_all_chr(tracts: pd.DataFrame,
                             short_cm: float = 1.0,
                             edge_support_cm: float = 2.0,
                             preview_n: int = 0,
                             sort_output: bool = True) -> pd.DataFrame:
    """
    Apply smooth_and_merge_chr() to each chromosome in a multi-chromosome
    haplotype tract dataframe and return the concatenated result.

    Parameters
    ----------
    tracts : DataFrame with columns including ['chm','sgpos','egpos','ancestry', ...]
    short_cm : float, threshold for a tract to be considered "short"
    edge_support_cm : float, minimum neighbor length at edges to allow a flip
    preview_n : int, lines to preview per-chrom changes (passed through; 0 = none)
    sort_output : bool, sort final result by ['chm','sgpos','egpos','id']

    Returns
    -------
    DataFrame with smoothed + merged tracts across all chromosomes.
    """
    if tracts.empty:
        return tracts.copy()

    chroms = tracts['chm'].unique()
    smoothed_list = []

    for chrom in chroms:
        print(f"Processing chromosome: {chrom}")
        tr_chrom = tracts.loc[tracts['chm'] == chrom].copy()
        tr_chrom_sm = smooth_and_merge_chr(
            tr_chrom,
            short_cm=short_cm,
            edge_support_cm=edge_support_cm,
            preview_n=preview_n
        )
        smoothed_list.append(tr_chrom_sm)

    out = pd.concat(smoothed_list, ignore_index=True) if smoothed_list else tracts.copy()
    if sort_output and {'chm','sgpos','egpos','id'}.issubset(out.columns):
        out = out.sort_values(['chm','sgpos','egpos','id']).reset_index(drop=True)
    return out


def main():
    parser = argparse.ArgumentParser(
        description="Smooth short noisy ancestry tracts per chromosome and merge consecutive same-ancestry runs."
    )
    parser.add_argument("input", help="Input tracts BED/TSV file (tab-separated)")
    parser.add_argument("output", help="Output BED/TSV file for smoothed+merged tracts")
    args = parser.parse_args()

    # read
    tracts = pd.read_csv(args.input, sep="\t")

    # defaults for quick testing
    minCM = 1
    edgeCM = 2.0

    # smooth + merge per chromosome
    smoothed = smooth_and_merge_all_chr(
        tracts, short_cm=minCM, edge_support_cm=edgeCM, preview_n=0
    )

    # write
    smoothed.to_csv(args.output, sep="\t", index=False)
    print(f"Done. Wrote {len(smoothed)} tracts to {args.output}")

if __name__ == "__main__":
    main()
