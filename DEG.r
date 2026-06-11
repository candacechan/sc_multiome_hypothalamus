suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(ggplot2)
})

params <- list(
  input_rds     = "path/to/seurat_object.RDS",
  celltype_col  = "celltype", 
  sample_col    = "orig.ident",
  out_dir       = "markers", 
  label         = "hypo_",

  assays        = c(RNA = "DEG", peaks = "DAR"),

  padj_cutoff   = 0.05,
  log2fc_cutoff = 1,

  n_top_heatmap = 1000,
  heatmap_w     = 40,
  heatmap_h     = 20
)

dir.create(params$out_dir, showWarnings = FALSE, recursive = TRUE)

get_counts <- function(obj, assay) {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, slot = "counts")
  )
}

run_pseudobulk <- function(seurat, assay, out_label, params) {
  bulk <- AggregateExpression(
    object        = seurat,
    assays        = assay,
    return.seurat = TRUE,
    group.by      = c(params$celltype_col, params$sample_col)
  )
  Idents(bulk) <- params$celltype_col

  de <- FindAllMarkers(object = bulk, test.use = "DESeq2", assay = assay)

  saveRDS(de, file.path(params$out_dir,
                        paste0("bulk_", out_label, "_", params$label, ".RDS")))
  write.csv(de, file.path(params$out_dir,
                          paste0("bulk_", out_label, "_", params$label, ".csv")),
            row.names = FALSE)

  sig <- de %>%
    filter(p_val_adj < params$padj_cutoff,
           abs(avg_log2FC) > params$log2fc_cutoff)
  write.csv(sig, file.path(params$out_dir,
                           paste0("bulk_", out_label, "_significant_",
                                  params$label, ".csv")),
            row.names = FALSE)

  top <- de %>%
    group_by(cluster) %>%
    filter(p_val_adj < params$padj_cutoff) %>%
    slice_max(n = params$n_top_heatmap, order_by = avg_log2FC)

  invisible(list(bulk = bulk, de = de, sig = sig))
}

seurat <- readRDS(params$input_rds)
stopifnot(params$celltype_col %in% colnames(seurat@meta.data),
          params$sample_col   %in% colnames(seurat@meta.data))
Idents(seurat) <- params$celltype_col

results <- lapply(names(params$assays), function(a)
  run_pseudobulk(seurat, assay = a, out_label = params$assays[[a]], params = params))
names(results) <- names(params$assays)
