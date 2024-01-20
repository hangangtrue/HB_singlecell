library(Seurat)

args <- unlist(strsplit(unlist(commandArgs()), " "))
temp.dir <- args[4]
### Mouse number
m <- as.numeric(args[5])
### Average cell number per mouse
ce <- as.numeric(args[6])
### Effect size, fold change
f <- as.numeric(args[7])
### Replicates
j <- as.numeric(args[8])
Hypo <- "Neuro_normal_chow.RData"
load(file.path(Hypo))
GoI <- c("Glp1r")
if(sum(table(combined$Sample_ID)>ce) <= 12 & m >= 12) m =10 # stop("Number of mice is smaller than 12. Stopping the script.")

### Available samples
id_aval <- names(table(combined$Sample_ID))[table(combined$Sample_ID)>(ce+100)]
combined <- combined %>%
  subset(Sample_ID %in% id_aval)

### Sample mouse from available samples
set.seed(j)
combined$status_p <- "HC"
s_m <- sample(id_aval, m)
combined <- combined %>%
  subset(Sample_ID %in% s_m)

### Assign m/2 to disease
s_m_d <- sample(s_m, m/2)
combined$status_p[combined$Sample_ID %in% s_m_d] <- "Disease"
combined$status_p  <- as.factor(combined$status_p)

### Sample cell numbers per mouse from Poisson
set.seed(ce+j)
n_ce <- rpois(m, lambda = ce)

### Subsample cells from each selected mouse
ind_1 <- c()
for(ii in 1:length(s_m)){
  ind <- which(combined$Sample_ID == s_m[ii])
  ind_1 <- c(ind_1, sample(ind, n_ce[ii]))
}
combined <- combined[, ind_1]

### Select genes expressed in at least 10% of cells and genes in GoI
pc <- 0.10
ind <- unique(c(which(rowMeans(combined@assays$RNA@counts != 0) > pc), which(rownames(combined) %in% GoI)))
combined <- combined[ind, ]

### Select 5% genes to add a fixed effect size
set.seed(m+ce+f+j)
n_g <- round(length(rownames(combined))*0.05)
s_g <- sample(rownames(combined), n_g)
s_g <- unique(c(s_g, GoI))
selected_tDE <- paste("selected_tDE_Mice_", m, "_Cells_", ce, "_FC_", f, "_Rep_",j,"_data.csv", sep = "")
write.csv(s_g, file = file.path("MAST",selected_tDE), row.names = FALSE)

### Add effect size (fold change) to the selected genes in the Disease group
combined[["RNA"]]@counts[rownames(combined) %in% s_g, combined$status_p == "Disease"] <- combined[["RNA"]]@counts[rownames(combined) %in% s_g, combined$status_p == "Disease"]*f

### Assign groups
combined$label <- combined$status_p
DefaultAssay(combined) <- "RNA"
meta_group = levels(combined$status_p)
ident.1 = meta_group[1]
ident.2 = meta_group[2]
Idents(combined) <- combined$C2_named

#### MAST
DEs.seurat <- function(i, combined, ident.1 = "disease", ident.2 = "nondisease", 
                       min.cells.feature = 0, min.cells.group = 0, logfc.threshold = 0, 
                       min.pct = 0, test.use = "MAST", verbose = TRUE, parallel = FALSE){
  
  combined$celltype.status <- paste(Idents(combined), combined$status_p, sep = "_")
  combined$celltype <- Idents(combined)
  Idents(combined) <- "celltype.status"
  
  clst <- i
  clst_status.1 <- paste(clst, ident.1, sep = "_")
  clst_status.2 <- paste(clst, ident.2, sep = "_")
  
  results <- FindMarkers(combined, ident.1 = clst_status.1, ident.2 = clst_status.2, 
                         min.cells.feature = min.cells.feature, min.cells.group = min.cells.group,
                         verbose = verbose, logfc.threshold = logfc.threshold, min.pct = min.pct, test.use = test.use)
  
  results <- results %>% mutate(cluster =i, gene = rownames(results))
  
  results
}

pct <- 0
min.cells.feature = 0
min.cells.group = 0
logfc.threshold = 0
test.use <- "MAST"
i <- "C2-1: Neurons"
DEs_file <- paste("Mice_", m, "_Cells_", ce, "_FC_", f, "_Rep_",j,"_data.csv", sep = "")

DE_genes <- DEs.seurat(i, combined, ident.1 = ident.1, ident.2 = ident.2, 
                       min.cells.feature = min.cells.feature, min.cells.group = min.cells.group, 
                       logfc.threshold = logfc.threshold, min.pct = pct,
                       test.use = test.use, verbose = TRUE, parallel = FALSE)

DE_genes <- DE_genes %>% 
  mutate(FC = 2^(avg_log2FC), FC = ifelse(FC<1, -1/FC, FC)) %>% 
  dplyr::select(p_val, FC, pct.1, pct.2, p_val_adj, cluster, gene) %>%
  arrange(p_val)

write.csv(DE_genes, file.path(temp.dir, DEs_file),row.names=F)
