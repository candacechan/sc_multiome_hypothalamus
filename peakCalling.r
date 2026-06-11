suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(BSgenome.Mmusculus.UCSC.mm10))

in_rds       <- "path/to/seurat_object.RDS"
group_column <- "celltype"
macs2_path   <- "/path/to/macs2"
out_dir      <- "MACS2_output"
genome       <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)

dir.create(out_dir, showWarnings = FALSE)

seurat <- readRDS(in_rds)
Idents(seurat) <- group_column

peaks <- CallPeaks(
  object           = seurat,
  assay            = "ATAC",
  group.by         = group_column,
  macs2.path       = macs2_path,
  outdir           = out_dir,
  fragment.tempdir = out_dir,
  cleanup          = FALSE
)

# Filter peaks to standard chromosomes and remove ENCODE blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

saveRDS(peaks, file = "peaks_celltype_specific.RDS")

DefaultAssay(seurat) <- "ATAC"
frags <- Fragments(seurat)

peak_counts <- FeatureMatrix(
  fragments = frags,
  features  = peaks,
  cells     = colnames(seurat)
)

seurat[["peaks"]] <- CreateChromatinAssay(
  counts     = peak_counts,
  sep        = c(":", "-"),
  genome     = genome,
  fragments  = frags,
  annotation = annotation
)

DefaultAssay(seurat) <- "peaks"
seurat <- RunTFIDF(seurat, method = 1)
seurat <- FindTopFeatures(seurat)

saveRDS(seurat, "celltype_peaks.RDS")
