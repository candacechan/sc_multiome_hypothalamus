library(Seurat)
library(Signac)
library(Pando)
library(igraph)
library(dplyr)
library(BSgenome)

work_dir <- "/path/to/analysis"
input_rds <- "seurat_object.rds"
label_col <- "broad_label"
genome_package <- "BSgenome.Hsapiens.UCSC.hg38"

setwd(work_dir)
dir.create("PANDO", showWarnings = FALSE)

genome <- getBSgenome(genome_package)
obj <- readRDS(input_rds)

labels <- as.character(obj@meta.data[[label_col]])
groups <- list(
    GLU = colnames(obj)[grepl("GLU", labels, ignore.case = TRUE) &
        !grepl("GABA", labels, ignore.case = TRUE)],
    GABA = colnames(obj)[grepl("GABA", labels, ignore.case = TRUE) &
        !grepl("GLU", labels, ignore.case = TRUE)]
)

data(motifs, package = "Pando")
data(motif2tf, package = "Pando")
motif2tf <- as.data.frame(motif2tf)

genes <- rownames(obj[["RNA"]])
gene_map <- setNames(genes, toupper(genes))
mapped_tf <- unname(gene_map[toupper(motif2tf$tf)])
motif2tf$tf[!is.na(mapped_tf)] <- mapped_tf[!is.na(mapped_tf)]

run_grn <- function(cells, name) {
    x <- obj[, cells]

    x <- initiate_grn(
        x,
        rna_assay = "RNA",
        peak_assay = "peaks",
        exclude_exons = FALSE
    )

    x <- find_motifs(
        x,
        pfm = motifs,
        motif_tfs = motif2tf,
        genome = genome
    )

    x <- infer_grn(
        x,
        peak_to_gene_method = "Signac",
        genes = genes[!grepl("^[0-9]", genes)],
        tf_cor = 0.1,
        parallel = TRUE
    )

    x <- find_modules(
        x,
        p_thresh = 0.1,
        nvar_thresh = 2,
        rsq_thresh = 0.05,
        min_genes_per_module = 1
    )

    x <- get_network_graph(x)
    saveRDS(x, file.path("PANDO", paste0(name, "_grn.rds")))
    x
}

grn <- lapply(names(groups), function(name) run_grn(groups[[name]], name))
names(grn) <- names(groups)

get_strength <- function(x, name) {
    net <- NetworkGraph(x)
    values <- strength(net, mode = "all", weights = abs(E(net)$weight))
    data.frame(
        gene = names(values),
        value = as.numeric(values),
        percentile = rank(values, ties.method = "average") / length(values)
    ) |>
        rename_with(~ paste0(name, "_", .x), -gene)
}

centrality <- full_join(
    get_strength(grn$GLU, "GLU"),
    get_strength(grn$GABA, "GABA"),
    by = "gene"
)

centrality[is.na(centrality)] <- 0
centrality$difference <- centrality$GLU_percentile - centrality$GABA_percentile
centrality <- centrality[order(abs(centrality$difference), decreasing = TRUE), ]

write.csv(
    centrality,
    file.path("PANDO", "GLU_vs_GABA_strength_centrality.csv"),
    row.names = FALSE
)
