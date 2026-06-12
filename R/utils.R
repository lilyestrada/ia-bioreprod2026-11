#' Download the yeast oscillating‑genes TSV (or any URL) safely
#'
#' @param url       URL of the file to download.
#' @param out_dir   Destination folder (default = "data").
#' @param overwrite Logical; if TRUE, forces re‑download.
#' @return          Full path to the downloaded file.
#' @export
# ---------------------------------------------------------------
# Simple downloader that wraps download.file()
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# download_yeast_file() – tiny wrapper around download.file()
# ---------------------------------------------------------------
download_yeast_file <- function(url,
                                out_dir = "data",
                                overwrite = FALSE) {
  # 1. make sure the destination folder exists
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 2. safe file name (strip possible query strings / fragments)
  file_name <- sub("[\\?#].*$", "", basename(url))
  if (nzchar(file_name) == FALSE) file_name <- "downloaded_file"
  
  dest <- file.path(out_dir, file_name)
  
  # 3. skip if we already have the file and overwrite = FALSE
  if (file.exists(dest) && !overwrite) {
    message("File already present – using cached copy: ", dest)
    return(normalizePath(dest, winslash = "/"))
  }
  
  # 4. download -------------------------------------------------
  message("Downloading ", url, "\n  → ", dest)
  
  # download.file() returns a status code (0 = success)
  # we wrap it in tryCatch so we can report a nice error if something goes wrong
  rc <- tryCatch(
    download.file(url, destfile = dest, mode = "wb", quiet = FALSE),
    error = function(e) e
  )
  
  # 5. error handling --------------------------------------------
  if (inherits(rc, "error") ||
      (!is.numeric(rc) && rc != 0) || rc != 0) {
    # remove a possibly half‑written file so a later retry starts clean
    if (file.exists(dest)) file.remove(dest)
    stop("Download failed (status = ", rc,
         "). Check the URL or your internet connection.")
  }
  
  message("Download finished.")
  normalizePath(dest, winslash = "/")
}


# Function to plot heatmaps
library(httr2)
library(jsonlite)
library(base64enc)

plot_heatmap_pdf <- function(outdir="output"){
  
  workdir <- getwd()
  
  analysis_dir <- "data"
  
  data_table <- file.path(analysis_dir, "oscillating-genes_1705_normalized-profiles.tsv")
  if (file.exists(data_table)) {
    message("data table\t", data_table)
  } else {
    message("First download the files with the script download_cell-cycle_files.R")
    stop("missing file\t", data_table)
  }
  
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
  #
  # Methods notes (from paper):
  #   - S. cerevisiae sampled every 5 min → 50 points = 0–245 min
  #   - Cells synchronised by alpha-factor arrest then released
  #   - Four periodicity algorithms ranked all genes; top 1600 kept
  #   - Lomb-Scargle p-value cutoff applied → final 1246 periodic genes
  #   - The TSV supplied already contains these 1246 pre-filtered genes
  #   - Expression depicted as z-score (SDs from per-gene mean)
  #   - Genes ordered on y-axis by peak expression time (Fig. 2 legend)
  #   - Colour scale: yellow = high (+1.5), black = 0, cyan = low (-1.5)
  # ============================================================
  
  # ---- 0. Required packages -----------------------------------
  # install.packages(c("gplots"))  # run once if needed
  library(gplots)
  
  # ---- 1. Load FPKM data --------------------------------------
  fpkm <- read.delim(
    data_table,
    row.names    = 1,
    check.names  = FALSE
  )
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  # Dimensions: 1705 genes × 50 samples (SRR accessions)
  
  # ---- 2. Assign time-point labels ----------------------------
  # 50 samples collected every 5 min starting at t = 0 min
  # (Methods: "every 5 minutes for S. cerevisiae")
  
  # ---- 3. Select the 1246 periodic genes ----------------------
  # The file contains 1705 candidate genes. The paper's three-step
  # filter (non-noisy AND top-1600 cumulative rank AND LS cutoff)
  # yields 1246. Because the exact per-gene LS p-values require the
  # supplementary S1 Table, the simplest reproducible approach is to
  # use the top 1246 rows as supplied (the file is already rank-ordered
  # by cumulative periodicity score; row order mirrors S1 Table col 2).
  #
  # If you have S1 Table, replace the line below with an explicit
  # gene list read from that file (column "LS_pass" == TRUE, top 1600).
  n_periodic <- 1246
  mat        <- as.matrix(fpkm[seq_len(n_periodic), ])
  
  # ---- 4. Z-score normalise each gene (row) -------------------
  # "Transcript levels depicted as z-score change relative to mean
  #  expression for each gene" (Fig. 2 legend)
  mat_z <- t(scale(t(mat)))   # scale() operates on columns → transpose twice
  
  # ---- 5. Order genes by peak expression time -----------------
  # "Genes ordered along y-axis by peak time of expression" (Fig. 2)
  peak_col   <- apply(mat_z, 1, which.max)
  mat_sorted <- mat_z[order(peak_col), ]
  
  # ---- 6. Colour palette: yellow → black → cyan ---------------
  # Matches the published figure exactly (high = yellow, 0 = black,
  # low = cyan)
  pal <- colorpanel(100, low = "cyan", mid = "black", high = "yellow")
  
  # ---- 7. Draw heatmap ----------------------------------------
  out_pdf <- file.path(outdir, "heatmap_yeast_cell-cycle.pdf")
  
  pdf(out_pdf, width = 7, height = 9)
  
  heatmap.2(
    mat_sorted,
    # --- layout ---
    Rowv         = FALSE,          # preserve peak-time row order
    Colv         = FALSE,          # preserve chronological column order
    dendrogram   = "none",
    # --- colour ---
    col          = pal,
    breaks       = seq(-1.5, 1.5, length.out = 101),
    # --- annotations ---
    trace        = "none",
    density.info = "none",
    key          = TRUE,
    keysize      = 1.2,
    key.title    = NA,
    key.xlab     = "",
    symkey       = TRUE,
    # --- labels ---
    labRow       = NA,             # 1246 gene names would be illegible
    labCol       = colnames(mat_sorted),
    cexCol       = 0.55,
    # --- margins & titles ---
    margins      = c(5, 5),
    xlab         = "time (minutes)",
    ylab         = paste0("Top Periodic Genes (", nrow(mat_sorted), ")"),
    main         = expression(italic("Saccharomyces cerevisiae"))
  )
  
  dev.off()
  
  message("Done — heatmap written to heatmap_yeast_cell-cycle.pdf")
  
}
