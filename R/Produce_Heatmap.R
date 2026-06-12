# ============================================================
# Atelier IA & Bioscripting 2026
# 12 juin 2026
# ============================================================

# Download data
DATA_URL <- "https://ifb-elixirfr.github.io/AI-for-scripting-bioanalysis/data/yeast-transcriptome-cell-cycle/oscillating-genes_normalized-profiles.tsv"
download_yeast_file(DATA_URL)

# Plot heatmap
plot_heatmap_pdf(outdir="output")
