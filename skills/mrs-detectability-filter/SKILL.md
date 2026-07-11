---
name: mrs-detectability-filter
description: >
  Given a perturbation or differential-expression dataset (e.g. Perturb-seq DE
  tables) and a set of candidate hits, decide which candidates could be
  validated in vivo with hyperpolarized 13C MRS. Filters hits to those that
  move a FLUX-DETERMINING node of a tracer-accessible reaction (not merely a
  pathway), predicts the change in the observable rate constant via a two-site
  exchange model, and flags candidates whose predicted change falls below the
  in vivo detectability floor. Use whenever someone asks "which of these
  metabolic hits could I actually measure with a hyperpolarized tracer / MRS /
  pyruvate imaging" or wants a pharmacodynamic-validation plan for a
  metabolic drug target.
license: MIT
---

# MRS Detectability Filter

## What this Skill does

Transcript fold-change != metabolic flux. Hyperpolarized 13C MRS is sensitive
to a *small number of rate-determining nodes*, not to whole pathways. This
Skill encodes that distinction so a genome-scale hit list can be reduced to the
candidates that are actually trackable in a living system.

## When to use it

- The user has candidate metabolic regulators (from a screen, a DE table, a
  literature list) and wants to know which are MRS-verifiable.
- The user wants a per-target pharmacodynamic-validation paragraph naming a
  specific tracer and an expected effect size.

## Inputs

1. A DE / effect-size table keyed by perturbation, with per-gene log2FC and a
   significance measure. (For the reference dataset: `GWCD4i.DE_stats.h5ad`,
   layers `log_fc`, `zscore`, `adj_p_value`.)
2. `config/gene_sets.yaml` — the Tier A flux nodes + tracer mapping. **Authored
   by the researcher.** The Skill refuses to invent this mapping.

## Pipeline

1. **Coverage check.** Intersect every gene set with tested genes; report
   found/requested per set. Abort a set below the coverage floor.
2. **Score** each perturbation against Tier B pathway sets and phenotype sets,
   masking the perturbed gene's own column.
3. **De-confound.** Regress each metabolic score on `n_downstream` (or the
   proliferation signature) and keep residuals — kills the ribosome/spliceosome
   artifact.
4. **Competitive null.** Empirical p vs. baseMean-matched random gene sets.
5. **Detectability filter.** Require >=1 Tier A node with |z|>=3, FDR<0.1, in a
   direction consistent with the Tier B shift. Record node, tracer, rate const.
6. **Forward model.** `src/forward_model.py` maps the node's log2FC to a
   predicted delta in the rate constant and applies the detectability threshold.
7. **Emit** a per-candidate package: effect sizes, guide/donor concordance,
   druggability tier, mechanism note, tracer + predicted delta + verdict.

## Outputs

- `candidates.tsv` — ranked, with detectability verdict per candidate.
- One markdown package per surviving candidate.
- A rejected-candidates list ("real biology, not PD-trackable") — keep this,
  it is the honest differentiator.

## Guardrails

- Never fabricate a tracer for a pathway that has no HP probe (e.g. one-carbon).
- Always report the transcript->flux `coupling` assumption; never present a
  predicted flux change as a measured one.
