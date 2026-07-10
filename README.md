# mrs-detectability-filter

**Find immunometabolic drug targets in a genome-scale T-cell CRISPRi screen, and
keep only the ones a hyperpolarized ¹³C MRS assay could actually verify in vivo.**

Built for *Built with Claude: Life Sciences* (Jul 7–13, 2026), Researcher track.

## The gap

Genome-scale perturbation screens produce metabolic and phenotypic hit lists
separately, with no built-in way to know which computational target could be
validated as a drug is developed against it — short of invasive biopsies or
endpoint assays. This project closes that gap for the metabolic axis.

## The idea

A perturbation is a valid hyperpolarized-MRS pharmacodynamic candidate **if and
only if** its knockdown moves a *flux-determining node* of a tracer-accessible
reaction — because HP-¹³C sees LDHA/LDHB + MCT1 (kPL) or the PDH complex (kPB),
not whole pathways. That converts "the metabolic program shifted" into a
physically necessary, testable filter. A two-site exchange forward model then
predicts the change in the observable rate constant and rejects any candidate
below the in vivo detectability floor.

## Repo layout

- `config/gene_sets.yaml` — Tier A flux nodes + tracer mapping, Tier B pathways,
  phenotype and proliferation sets. **Authored by the researcher.**
- `src/forward_model.py` — HP-¹³C two-site exchange model + detectability verdict.
- `skills/mrs-detectability-filter/SKILL.md` — the reusable Claude Science Skill.

## Run

```bash
pip install scipy numpy
python src/forward_model.py    # toy LDHA-knockdown demo
```

## Data

Primary Human CD4+ T Cell Perturb-seq (CRISPRi), Zhu, Dann, Yan, Reyes Retana,
Goto, Guitche, Petersen, Ota, Pritchard, Marson.
bioRxiv 2025.12.23.696273. Dataset: MIT License.
Analysis reference: https://github.com/emdann/GWT_perturbseq_analysis_2025

Files used: `GWCD4i.DE_stats.h5ad`, `GWCD4i.DE_stats.by_guide.h5mu`,
`DE_stats.suppl_table.csv`.

## Method provenance

Hyperpolarized ¹³C MRS of activated human T-cell metabolic reprogramming:
Can et al., *Scientific Reports*, 2020 (cited as prior methodological work,
not as project content).

## License

MIT — see `LICENSE`.
