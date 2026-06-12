# AI-bioscripting-figure

  # ============================================================
  # Yeast Cell Cycle Heatmap — Saccharomyces cerevisiae
  # Reproduces Fig. 2A from:
  #   Kelliher et al. (2016) PLoS Genetics 12(12): e1006453
  #   doi:10.1371/journal.pgen.1006453
  #
  # Input : oscillating-genes_1705_normalized-profiles.tsv
  #         (1705 periodic-gene candidates × 50 time points;
  #          CuffNorm-normalised FPKM values; GEO: GSE80474)
  #
  # Output: heatmap_yeast_cell-cycle.pdf

Use case:
# Download data
DATA_URL <- "https://ifb-elixirfr.github.io/AI-for-scripting-bioanalysis/data/yeast-transcriptome-cell-cycle/oscillating-genes_normalized-profiles.tsv"
download_yeast_file(DATA_URL)

# Plot heatmap
plot_heatmap_pdf(outdir="output")
