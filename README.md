# AI-bioscripting-figure

## Yeast Cell Cycle Heatmap (*Saccharomyces cerevisiae*)

Reproduces **Figure 2A** from:

> Kelliher et al. (2016). *Improved Deconvolution of Cell Cycle Population Data Reveals Gene Expression Patterns in Saccharomyces cerevisiae*. PLoS Genetics, 12(12): e1006453.  
> DOI: 10.1371/journal.pgen.1006453

### Input

`oscillating-genes_normalized-profiles.tsv`

- 1,705 periodic gene candidates
- 50 cell-cycle time points
- CuffNorm-normalized FPKM expression values
- Source dataset: GEO accession **GSE80474**

### Output

`heatmap_yeast_cell-cycle.pdf`

Heatmap showing normalized expression profiles of oscillating genes throughout the yeast cell cycle.

## Use Case

```r
# Download data
DATA_URL <- "https://ifb-elixirfr.github.io/AI-for-scripting-bioanalysis/data/yeast-transcriptome-cell-cycle/oscillating-genes_normalized-profiles.tsv"
download_yeast_file(DATA_URL)

# Plot heatmap
plot_heatmap_pdf(outdir = "output")
```

The generated PDF heatmap is saved in the `output/` directory.
