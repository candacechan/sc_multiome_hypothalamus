library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(S4Vectors)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(presto)
library(glue)
library(scales)
library(enrichR)
library(GenomeInfoDb)

seurat <- readRDS("subset_diet.RDS")

astrocyte_labels <- unique(seurat$C2_named_label[
  grepl("Astrocyte", seurat$C2_named_label, ignore.case = TRUE)
])

seurat_astros <- subset(seurat, subset = C2_named_label %in% astrocyte_labels)
Idents(seurat_astros) <- seurat_astros$C2_named_label

top_markers <- FindAllMarkers(
  object = seurat_astros,
  test.use = "wilcox",
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.05,
  logfc.threshold = 0.25
)

plot_enrich <- function(cluster){
    topDEG <- top_markers %>%
        filter(cluster == !!cluster) %>%
        filter(p_val_adj < 0.05) %>% 
        filter(avg_log2FC > 1) %>% 
        slice_max(n = 1000, order_by = avg_log2FC)
      
    neg.markers.list <- unique(topDEG$gene)
    enrich.database <- "GO_Molecular_Function_2023"
    num.pathway <- 5
    fontsize <- 8
    
    # Run enrichR
    neg.er <- enrichR::enrichr(genes = neg.markers.list, databases = enrich.database)
    neg.er <- neg.er[[enrich.database]] # Cleaner way to access the list
    
    if (nrow(neg.er) == 0) {
        message(paste("No enrichment found for:", cluster))
        return(NULL)
    }

    neg.er$log10pval <- -log10(neg.er$P.value)
    neg.er <- neg.er[1:min(num.pathway, nrow(neg.er)), ]
    neg.er$term <- factor(neg.er$Term, levels = rev(neg.er$Term)) # Rev for better coord_flip order

    p2 <- ggplot(data = neg.er, aes(x = term, y = log10pval)) + 
        geom_bar(stat = "identity", fill = "darkolivegreen3") + 
        coord_flip() + 
        xlab("Pathway") + 
        ylab("-log10(pval)") + 
        ggtitle(paste(cluster, "GO MF")) + 
        theme_classic() + 
        theme(text = element_text(size = fontsize))
    
    print(p2)
}

lapply(astrocyte_labels, plot_enrich)
