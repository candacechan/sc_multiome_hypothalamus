### Code used for analysis of multiome (RNA+ATAC) for mouse dataset
### C. Chan

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(S4Vectors)
library(ggplot2)
library(patchwok)
library(GenomicRanges)
library(future)
library(dplyr)
library(harmony)
library(presto)

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"
genome <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)

counts_h1 <- Read10X_h5("M1/outs/filtered_feature_bc_matrix.h5")
counts_h2 <- Read10X_h5("M2/outs/filtered_feature_bc_matrix.h5")
counts_h3 <- Read10X_h5("M3/outs/filtered_feature_bc_matrix.h5")
counts_h4 <- Read10X_h5("M4/outs/filtered_feature_bc_matrix.h5")
counts_h5 <- Read10X_h5("M5/outs/filtered_feature_bc_matrix.h5")
counts_h6 <- Read10X_h5("M6/outs/filtered_feature_bc_matrix.h5")

seurat_h1 <- CreateSeuratObject(counts = counts_h1$`Gene Expression`,
										assay = "RNA",
										project = "m1")
seurat_h2 <- CreateSeuratObject(counts = counts_h2$`Gene Expression`,
										assay = "RNA",
										project = "m2")
seurat_h3 <- CreateSeuratObject(counts = counts_h3$`Gene Expression`,
										assay = "RNA",
										project = "m3")
seurat_h4 <- CreateSeuratObject(counts = counts_h4$`Gene Expression`,
										assay = "RNA",
										project = "m4")
seurat_h5 <- CreateSeuratObject(counts = counts_h5$`Gene Expression`,
										assay = "RNA",
										project = "m5")
seurat_h6 <- CreateSeuratObject(counts = counts_h6$`Gene Expression`,
										assay = "RNA",
										project = "m6")

seurat_h1[['ATAC']] <- CreateChromatinAssay(counts = counts_h1$`Peaks`,
											 annotation = annotation,
											 fragments = "M1/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome)
seurat_h2[['ATAC']] <- CreateChromatinAssay(counts = counts_h2$`Peaks`,
											 annotation = annotation,
											 fragments = "M2/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome)
seurat_h3[['ATAC']] <- CreateChromatinAssay(counts = counts_h3$`Peaks`,
											 annotation = annotation,
											 fragments = "M3/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome)
seurat_h4[['ATAC']] <- CreateChromatinAssay(counts = counts_h4$`Peaks`,
											 annotation = annotation,
											 fragments = "M4/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome)
seurat_h5[['ATAC']] <- CreateChromatinAssay(counts = counts_h5$`Peaks`,
											 annotation = annotation,
											 fragments = "M5/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome)											 
seurat_h6[['ATAC']] <- CreateChromatinAssay(counts = counts_h6$`Peaks`,
											 annotation = annotation,
											 fragments = "M6/outs/atac_fragments.tsv.gz",
											 sep = c(":", "-"),
											 genome = genome) 

seurat <- merge(seurat_h1, c(seurat_h2, seurat_h3, seurat_h4, seurat_h5, seurat_h6))

seurat <- PercentageFeatureSet(seurat, pattern = "^mt-", col.name = "percent.mt", assay = "RNA")
seurat <- NucleosomeSignal(seurat, assay = "ATAC")
seurat <- TSSEnrichment(seurat, assay = "ATAC")

seurat <- subset(seurat,
  subset = nFeature_RNA > 1000 &
	nFeature_RNA < 7500 &
	percent.mt < 30 &
	nFeature_ATAC > 1000 &
	nFeature_ATAC < 30000 &
	TSS.enrichment > 1 &
	nucleosome_signal < 2
)

peaks <- CallPeaks(seurat,
				   assay="ATAC",
				   macs2.path="~/miniconda3/envs/macs/bin/macs2")
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

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
  RunHarmony(group.by.vars = "orig.ident", dims.use = 1:20, max.iter.harmony = 50, reduction.save = "harmony_rna") %>%
  RunUMAP(reduction = "harmony_rna", reduction.key = "UMAPRNA_", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony_rna", dims = 1:20) %>% 
  FindClusters(resolution = 0.6)

DefaultAssay(seurat) <- "ATAC"
seurat <- FindTopFeatures(seurat, min.cutoff=3)
seurat <- RunTFIDF(seurat, method=1)
seurat <- RunSVD(seurat, reduction.name="lsi", n=50)

seurat <- RunHarmony(seurat,
					 group.by.vars = "orig.ident",
					 reduction = "lsi",
					 dims.use = 2:30,
					 max.iter.harmony = 50,
					 reduction.save = "harmony_atac")
seurat <- RunUMAP(seurat,
				  reduction = "harmony_atac",
				  dims = 1:ncol(Embeddings(seurat,"harmony_atac")),
				  reduction.name = "umap_harmony_atac",
				  reduction.key = "UMAPHARMONYATAC_")

resolutions <- c(0.2, 0.4, 0.6)

seurat <- FindMultiModalNeighbors(seurat,
								  reduction.list = list("harmony_rna", "harmony_atac"),
								  dims.list = list(1:ncol(Embeddings(seurat,"harmony_rna")), 1:ncol(Embeddings(seurat,"harmony_atac"))),
								  modality.weight.name = c("RNA.weight","ATAC.weight"),
								  verbose = TRUE)

seurat <- RunUMAP(seurat, nn.name = "weighted.nn", 
				reduction.name = "wnn.umap", 
				reduction.key = "wnnUMAP_")
seurat <- FindClusters(seurat, graph.name = "wsnn", resolution = resolutions)
