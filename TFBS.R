library(Seurat)
library(Signac)
library(dplyr)
library(TFBSTools)
library(JASPAR2020)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(pheatmap)
  
input_rds <- "seurat_object.RDS"
label_col <- "broad_label"
out_dir <- "TF"
top_n <- 5

dir.create(out_dir, showWarnings = FALSE)

obj <- readRDS(input_rds)

pfm <- getMatrixSet(
    JASPAR2020,
    opts = list(species = 9606, all_versions = FALSE)
)

DefaultAssay(obj) <- "peaks"
obj <- AddMotifs(
    object = obj,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
)
obj <- RunChromVAR(
    object = obj,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(obj) <- "chromvar"
Idents(obj) <- obj@meta.data[[label_col]]

motif_markers <- FindAllMarkers(
    object = obj,
    test.use = "wilcox",
    only.pos = FALSE,
    min.pct = 0.1,
    mean.fxn = rowMeans,
    fc.name = "avg_diff"
)

saveRDS(motif_markers, file.path(out_dir, "motif_markers.rds"))
write.csv(
    motif_markers,
    file.path(out_dir, "motif_markers.csv"),
    row.names = FALSE
)

top_motifs <- motif_markers %>%
    filter(p_val_adj < 0.05, avg_diff > 0) %>%
    group_by(cluster) %>%
    slice_max(avg_diff, n = top_n, with_ties = FALSE) %>%
    ungroup()

features <- unique(top_motifs$gene)

avg_z <- AverageExpression(
    obj,
    assays = "chromvar",
    features = features,
    group.by = label_col,
    slot = "data"
)[["chromvar"]]

write.csv(
    avg_z,
    file.path(out_dir, "top_motif_average_chromvar_zscores.csv")
)

pheatmap(
    avg_z,
    cluster_cols = FALSE,
    border_color = NA,
    filename = file.path(out_dir, "motif_heatmap.pdf"),
    width = 5,
    height = 6
)
