library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(S4Vectors)
library(GenomicRanges)
library(JASPAR2020)
library(TFBSTools)
library(patchwork)

### LinkPeaks
DefaultAssay(seurat) <- "ATAC"
seurat <- RegionStats(seurat, genome = BSgenome.Mmusculus.UCSC.mm10)
seurat <- LinkPeaks(
	object = seurat,
	peak.assay = "ATAC",
	expression.assay = "RNA",
	genes.use = genes_of_interest,
	min.cells = 3,
	pvalue_cutoff = 0.05,
	distance = 1e+06
)

### DEG, DAR
DefaultAssay(seurat) <- "RNA"
markers<- FindAllMarkers(
						object = seurat, 
						only.pos = TRUE, 
						min.pct = 0.1,
						verbose=TRUE
						)
DefaultAssay(seurat) <- "ATAC"
da_peaks <- FindAllMarkers(
						object = seurat, 
						only.pos = TRUE, 
						min.pct = 0.05,
						latent.vars = 'nCount_ATAC')

### TFBS
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

seurat <- RunChromVAR(seurat, genome = BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(seurat) <- "chromvar"

TF_motifs <- FindAllMarkers(
							object = seurat, 
							only.pos = TRUE, 
							min.pct = 0.1,
							mean.fxn = rowMeans,
							fc.name = "avg_diff"
							)

### Sex-differential regions and DEG
genesSO <- subset(x = seurat, orig.ident == "h1" |orig.ident == "h2"| orig.ident == "h3")
seurat$sex <- ifelse(colnames(seurat) %in% colnames(genesSO), "M", "F")
seurat$annotated_mf <- paste(Idents(seurat), seurat$sex, sep="_")
Idents(seurat) <- "annotated_mf"

DefaultAssay(seurat) <- "RNA"
mf_deg <- FindAllMarkers(
						object = seurat, 
						only.pos = TRUE, 
						min.pct = 0.1,
						verbose=TRUE
						)

DefaultAssay(seurat) <- "ATAC"
mf_peaks <- FindAllMarkers(
						object = seurat, 
						only.pos = TRUE, 
						min.pct = 0.05,
						test.use = 'LR',
						verbose=TRUE
						latent.vars = 'nCount_ATAC')
