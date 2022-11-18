library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(S4Vectors)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
library(dplyr)
library(harmony)
library(presto)
library(SoupX)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
genome <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

counts_h1 <- Read10X_h5("Hypodata-0821/human_hypo_1/outs/filtered_feature_bc_matrix.h5")
counts_h2 <- Read10X_h5("Hypodata-0821/human_hypo_2/outs/filtered_feature_bc_matrix.h5")
counts_h3 <- Read10X_h5("Hypodata-0821/human_hypo_3/outs/filtered_feature_bc_matrix.h5")
counts_h4 <- Read10X_h5("Hypodata-0821/human_hypo_4/outs/filtered_feature_bc_matrix.h5")
counts_h5 <- Read10X_h5("Hypodata-0821/human_hypo_5/outs/filtered_feature_bc_matrix.h5")
counts_h6 <- Read10X_h5("Hypodata-0821/human_hypo_6/outs/filtered_feature_bc_matrix.h5")
counts_h7 <- Read10X_h5("Hypodata-091422_cellrangerOutput/M1/outs/filtered_feature_bc_matrix.h5")
counts_h8 <- Read10X_h5("Hypodata-091422_cellrangerOutput/M2/outs/filtered_feature_bc_matrix.h5")
counts_h9 <- Read10X_h5("Hypodata-091422_cellrangerOutput/M3/outs/filtered_feature_bc_matrix.h5")
counts_h10 <- Read10X_h5("Hypodata-091422_cellrangerOutput/M4/outs/filtered_feature_bc_matrix.h5")
counts_h11 <- Read10X_h5("Hypodata-091422_cellrangerOutput/M5/outs/filtered_feature_bc_matrix.h5")
counts_h12 <- Read10X_h5("Hypodata-091422_cellrangerOutput/M6/outs/filtered_feature_bc_matrix.h5")

sc_h1 = load10X("Hypodata-0821/human_hypo_1/outs")
sc_h2 = load10X("Hypodata-0821/human_hypo_2/outs")
sc_h3 = load10X("Hypodata-0821/human_hypo_3/outs")
sc_h4 = load10X("Hypodata-0821/human_hypo_4/outs")
sc_h5 = load10X("Hypodata-0821/human_hypo_5/outs")
sc_h6 = load10X("Hypodata-0821/human_hypo_6/outs")

sc_h7 = load10X("Hypodata-091422_cellrangerOutput/M1/outs/")
sc_h8 = load10X("Hypodata-091422_cellrangerOutput/M2/outs/")
sc_h9 = load10X("Hypodata-091422_cellrangerOutput/M3/outs/")
sc_h10 = load10X("Hypodata-091422_cellrangerOutput/M4/outs/")
sc_h11 = load10X("Hypodata-091422_cellrangerOutput/M5/outs/")
sc_h12 = load10X("Hypodata-091422_cellrangerOutput/M6/outs/")

# Estimate rho
sc_h1 = autoEstCont(sc_h1)
sc_h2 = autoEstCont(sc_h2)
sc_h3 = autoEstCont(sc_h3)
sc_h4 = autoEstCont(sc_h4)
sc_h5 = autoEstCont(sc_h5)
sc_h6 = autoEstCont(sc_h6)
sc_h7 = autoEstCont(sc_h7)
sc_h8 = autoEstCont(sc_h8)
sc_h9 = autoEstCont(sc_h9)
sc_h10 = autoEstCont(sc_h10)
sc_h11 = autoEstCont(sc_h11)
sc_h12 = autoEstCont(sc_h12)

# Clean the data
out_h1 = adjustCounts(sc_h1)
out_h2 = adjustCounts(sc_h2)
out_h3 = adjustCounts(sc_h3)
out_h4 = adjustCounts(sc_h4)
out_h5 = adjustCounts(sc_h5)
out_h6 = adjustCounts(sc_h6)
out_h7 = adjustCounts(sc_h7)
out_h8 = adjustCounts(sc_h8)
out_h9 = adjustCounts(sc_h9)
out_h10 = adjustCounts(sc_h10)
out_h11 = adjustCounts(sc_h11)
out_h12 = adjustCounts(sc_h12)

seurat_h1 <- CreateSeuratObject(counts = out_h1,
								assay = "RNA",
								project = "h1")
seurat_h2 <- CreateSeuratObject(counts = out_h2,
								assay = "RNA",
								project = "h2")
seurat_h3 <- CreateSeuratObject(counts = out_h3,
								assay = "RNA",
								project = "h3")
seurat_h4 <- CreateSeuratObject(counts = out_h4,
								assay = "RNA",
								project = "h4")
seurat_h5 <- CreateSeuratObject(counts = out_h5,
								assay = "RNA",
								project = "h5")
seurat_h6 <- CreateSeuratObject(counts = out_h6,
								assay = "RNA",
								project = "h6")
seurat_h7 <- CreateSeuratObject(counts = out_h7,
								assay = "RNA",
								project = "h7")
seurat_h8 <- CreateSeuratObject(counts = out_h8,
								assay = "RNA",
								project = "h8")
seurat_h9 <- CreateSeuratObject(counts = out_h9,
								assay = "RNA",
								project = "h9")
seurat_h10 <- CreateSeuratObject(counts = out_h10,
								assay = "RNA",
								project = "h10")
seurat_h11 <- CreateSeuratObject(counts = out_h11,
								assay = "RNA",
								project = "h11")
seurat_h12 <- CreateSeuratObject(counts = out_h12,
								assay = "RNA",
								project = "h12")

seurat_h1[['ATAC']] <- CreateChromatinAssay(counts = counts_h1$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-0821/human_hypo_1/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h2[['ATAC']] <- CreateChromatinAssay(counts = counts_h2$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-0821/human_hypo_2/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h3[['ATAC']] <- CreateChromatinAssay(counts = counts_h3$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-0821/human_hypo_3/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h4[['ATAC']] <- CreateChromatinAssay(counts = counts_h4$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-0821/human_hypo_4/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h5[['ATAC']] <- CreateChromatinAssay(counts = counts_h5$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-0821/human_hypo_5/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h6[['ATAC']] <- CreateChromatinAssay(counts = counts_h6$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-0821/human_hypo_6/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)							 
seurat_h7[['ATAC']] <- CreateChromatinAssay(counts = counts_h7$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-091422_cellrangerOutput/M1/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h8[['ATAC']] <- CreateChromatinAssay(counts = counts_h8$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-091422_cellrangerOutput/M2/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h9[['ATAC']] <- CreateChromatinAssay(counts = counts_h9$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-091422_cellrangerOutput/M3/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h10[['ATAC']] <- CreateChromatinAssay(counts = counts_h10$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-091422_cellrangerOutput/M4/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h11[['ATAC']] <- CreateChromatinAssay(counts = counts_h11$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-091422_cellrangerOutput/M5/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)
seurat_h12[['ATAC']] <- CreateChromatinAssay(counts = counts_h12$`Peaks`,
											 annotation = annotation,
											 fragments = "Hypodata-091422_cellrangerOutput/M6/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome,
											)	

seurat <- merge(seurat_h1, c(seurat_h2, seurat_h3, seurat_h4, seurat_h5, seurat_h6, seurat_h7, seurat_h8, seurat_h9, seurat_h10, seurat_h11, seurat_h12))

# # # # # # QC for RNA and ATAC
seurat <- PercentageFeatureSet(seurat, pattern = "^MT-", col.name = "percent.mt", assay = "RNA")
seurat <- NucleosomeSignal(seurat, assay = "ATAC")
seurat <- TSSEnrichment(seurat, assay = "ATAC")

seurat$sampleID <- seurat$orig.ident
Idents(object = seurat) <- seurat@meta.data$`sampleID`
seurat <- RenameIdents(seurat, `h7` = "h1", `h8` = "h2", `h9` = "h3", `h10` = "h4", `h11` = "h5", `h12` = "h6")
seurat$sampleID <- Idents(seurat)

seurat <- subset(seurat,
  subset = nFeature_RNA > 1000 &
	nFeature_RNA < 7500 &
	percent.mt < 30 &
	nFeature_ATAC > 1000 &
	nFeature_ATAC < 20000 &
	TSS.enrichment > 1 &
	nucleosome_signal < 2
)

peaks <- CallPeaks(seurat,
				   assay="ATAC",
				   macs2.path="~/miniconda3/envs/macs/bin/macs2")
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

counts_atac <- FeatureMatrix(fragments = seurat@assays$ATAC@fragments,
							features = peaks,
							cells = colnames(seurat)
							)				   

seurat[['ATAC']] <- CreateChromatinAssay(counts = counts_atac,
										 fragments = seurat@assays$ATAC@fragments,
										 annotation = seurat@assays$ATAC@annotation,
										 sep = c(":","-"),
										 genome = genome)

DefaultAssay(seurat) <- "RNA"

seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData(vars.to.regress = c("nFeature_RNA", "percent.mt")) %>%
  RunPCA(npcs = 50) %>%
  RunHarmony(group.by.vars = c("orig.ident", "sampleID"), dims.use = 1:20, max.iter.harmony = 50, reduction.save = "harmony_rna") %>%
  RunUMAP(reduction = "harmony_rna", reduction.key = "UMAPRNA_", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony_rna", dims = 1:20) %>% 
  FindClusters(resolution = c(0.4))

DefaultAssay(seurat) <- "ATAC"
seurat <- FindTopFeatures(seurat, min.cutoff=3)
seurat <- RunTFIDF(seurat, method=1)
seurat <- RunSVD(seurat, reduction.name="lsi", n=50)

seurat <- RunHarmony(seurat,
					 group.by.vars = c("orig.ident", "sampleID"),
					 reduction = "lsi",
					 dims.use = 2:30,
					 max.iter.harmony = 50,
					 reduction.save = "harmony_atac")
seurat <- RunUMAP(seurat,
				  reduction = "harmony_atac",
				  dims = 1:ncol(Embeddings(seurat,"harmony_atac")),
				  reduction.name = "umap_harmony_atac",
				  reduction.key = "UMAPHARMONYATAC_")

DefaultAssay(seurat) <- "RNA"

resolutions <- c(0.2, 0.4, 0.6, 0.8)

seurat <- FindMultiModalNeighbors(seurat,
								reduction.list = list("harmony_rna", "harmony_atac"),
								dims.list = list(
								1:ncol(Embeddings(seurat,"harmony_rna")),
								1:ncol(Embeddings(seurat,"harmony_atac"))),
								modality.weight.name = c("RNA.weight","ATAC.weight")
								)

seurat <- RunUMAP(seurat, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seurat <- FindClusters(seurat, graph.name = "wsnn", resolution = resolutions)
