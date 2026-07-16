library(Seurat)
library(dplyr)
library(tidyr)
library(Matrix)

input_rds <- "path/to/seurat_object.rds"
output_dir <- "sex_differential_analysis"

celltype_col <- "broad_label"
sample_col <- "sampleID"
sex_col <- "sex"

assays_to_test <- c("RNA", "peaks")
min_cells_per_pseudobulk <- 10
min_samples_per_sex <- 3

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(input_rds)

required_cols <- c(celltype_col, sample_col, sex_col)
missing_cols <- setdiff(required_cols, colnames(obj@meta.data))
if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

obj$analysis_celltype <- as.character(obj@meta.data[[celltype_col]])
obj$analysis_sample <- as.character(obj@meta.data[[sample_col]])
obj$analysis_sex <- as.character(obj@meta.data[[sex_col]])

obj$analysis_sex[obj$analysis_sex %in% c("Male", "male")] <- "M"
obj$analysis_sex[obj$analysis_sex %in% c("Female", "female")] <- "F"

sample_sex <- obj@meta.data %>%
    distinct(analysis_sample, analysis_sex)

if (anyDuplicated(sample_sex$analysis_sample)) {
    stop("Each biological sample must map to exactly one sex.")
}

cell_counts <- obj@meta.data %>%
    count(
        analysis_celltype,
        analysis_sample,
        analysis_sex,
        name = "n_cells"
    )

eligible_samples <- cell_counts %>%
    filter(n_cells >= min_cells_per_pseudobulk)

eligible_celltypes <- eligible_samples %>%
    count(analysis_celltype, analysis_sex, name = "n_samples") %>%
    pivot_wider(
        names_from = analysis_sex,
        values_from = n_samples,
        values_fill = 0
    ) %>%
    filter(M >= min_samples_per_sex, F >= min_samples_per_sex) %>%
    pull(analysis_celltype)

if (length(eligible_celltypes) == 0) {
    stop("No cell types meet the representation thresholds.")
}

obj <- subset(obj, subset = analysis_celltype %in% eligible_celltypes)

build_pseudobulk <- function(object, assay_name) {
    bulk <- AggregateExpression(
        object = object,
        assays = assay_name,
        group.by = c("analysis_celltype", "analysis_sample"),
        slot = "counts",
        return.seurat = TRUE
    )

    DefaultAssay(bulk) <- assay_name

    celltypes <- sort(unique(object$analysis_celltype), decreasing = TRUE)
    samples <- unique(object$analysis_sample)

    parsed_celltype <- rep(NA_character_, ncol(bulk))
    parsed_sample <- rep(NA_character_, ncol(bulk))

    for (i in seq_along(colnames(bulk))) {
        column_name <- colnames(bulk)[i]

        celltype_hit <- celltypes[
            startsWith(column_name, paste0(celltypes, "_"))
        ]
        sample_hit <- samples[
            endsWith(column_name, paste0("_", samples))
        ]

        if (length(celltype_hit) > 0) {
            parsed_celltype[i] <- celltype_hit[1]
        }
        if (length(sample_hit) > 0) {
            parsed_sample[i] <- sample_hit[1]
        }
    }

    if (anyNA(parsed_celltype) || anyNA(parsed_sample)) {
        stop("Unable to parse one or more pseudobulk column names.")
    }

    bulk$analysis_celltype <- parsed_celltype
    bulk$analysis_sample <- parsed_sample
    bulk$sex <- factor(
        sample_sex$analysis_sex[
            match(bulk$analysis_sample, sample_sex$analysis_sample)
        ],
        levels = c("F", "M")
    )

    eligible_keys <- eligible_samples %>%
        filter(analysis_celltype %in% eligible_celltypes) %>%
        transmute(
            key = paste(
                analysis_celltype,
                analysis_sample,
                sep = "___"
            )
        ) %>%
        pull(key)

    bulk_keys <- paste(
        bulk$analysis_celltype,
        bulk$analysis_sample,
        sep = "___"
    )

    bulk[, bulk_keys %in% eligible_keys]
}

run_sex_analysis <- function(object, assay_name) {
    if (!assay_name %in% Assays(object)) {
        warning("Skipping missing assay: ", assay_name)
        return(NULL)
    }

    bulk <- build_pseudobulk(object, assay_name)

    if (assay_name == "peaks") {
        counts <- GetAssayData(
            bulk,
            assay = assay_name,
            slot = "counts"
        )

        keep <- Matrix::rowSums(counts) >= 10 &
            Matrix::rowSums(counts > 0) >= 3

        bulk <- subset(
            bulk,
            features = rownames(counts)[keep]
        )
    }

    results <- list()

    for (celltype in sort(unique(bulk$analysis_celltype))) {
        celltype_bulk <- bulk[, bulk$analysis_celltype == celltype]
        Idents(celltype_bulk) <- "sex"

        n_male <- sum(celltype_bulk$sex == "M")
        n_female <- sum(celltype_bulk$sex == "F")

        if (
            n_male < min_samples_per_sex ||
            n_female < min_samples_per_sex
        ) {
            next
        }

        result <- tryCatch(
            FindMarkers(
                object = celltype_bulk,
                assay = assay_name,
                slot = "counts",
                ident.1 = "M",
                ident.2 = "F",
                test.use = "DESeq2",
                min.pct = 0,
                logfc.threshold = 0,
                min.cells.group = 1,
                verbose = FALSE
            ),
            error = function(e) NULL
        )

        if (is.null(result) || nrow(result) == 0) {
            next
        }

        result$feature <- rownames(result)
        result$cell_type <- celltype
        result$n_male_samples <- n_male
        result$n_female_samples <- n_female
        result$contrast <- "M_vs_F"

        results[[celltype]] <- result
    }

    if (length(results) == 0) {
        return(NULL)
    }

    bind_rows(results)
}

for (assay_name in assays_to_test) {
    result <- run_sex_analysis(obj, assay_name)

    if (is.null(result)) {
        next
    }

    prefix <- if (assay_name == "RNA") {
        "sex_DEG"
    } else {
        "sex_DAR"
    }

    saveRDS(
        result,
        file.path(output_dir, paste0(prefix, ".rds"))
    )

    write.csv(
        result,
        file.path(output_dir, paste0(prefix, ".csv")),
        row.names = FALSE
    )
}
