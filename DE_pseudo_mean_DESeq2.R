library(hierarchicell)
library(nebula)
library(Libra)
library(pbmcapply)

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
write.csv(s_g, file = file.path("pseudo_mean_DESeq2",selected_tDE), row.names = FALSE)

### Add effect size (fold change) to the selected genes in the Disease group
combined[["RNA"]]@counts[rownames(combined) %in% s_g, combined$status_p == "Disease"] <- combined[["RNA"]]@counts[rownames(combined) %in% s_g, combined$status_p == "Disease"]*f

### Assign groups
combined$label <- combined$status_p
DefaultAssay(combined) <- "RNA"
meta_group = levels(combined$status_p)
ident.1 = meta_group[1]
ident.2 = meta_group[2]
Idents(combined) <- combined$C2_named

##### Option 2 in the hierarchicell package (Pseudo-bulk mean)
i <- "C2-1: Neurons"
test.use <- "DESeq2_mean"
DEs_file <- paste("Mice_", m, "_Cells_", ce, "_FC_", f, "_Rep_",j,"_data.csv", sep = "")

combined$Sample_ID <- gsub("_", "-", combined$Sample_ID)
meta <- data.frame(wellKey = colnames(combined), DonorID = paste(combined$status_p, combined$Sample_ID, sep = "_"), Status = combined$status_p)
genecounts_degs <- t(combined@assays$RNA@counts)
genecounts_degs <- cbind(meta[,1:2],genecounts_degs)
computemeans_degs <- function(x){tapply(x,genecounts_degs[,2],mean)}
cellmeans_degs <- sapply(genecounts_degs[,c(-1,-2)],computemeans_degs)
rownames(cellmeans_degs) <- paste0(rownames(cellmeans_degs),"_Mean")
coldata_degs <- as.data.frame(cbind(rownames(cellmeans_degs),rownames(cellmeans_degs)))
colnames(coldata_degs) <- c("SampleID","ToSep")
coldata_degs <- tidyr::separate(coldata_degs,ToSep,c("Status", "Donor_Number", "Mean"), sep="_")
rownames(coldata_degs) <- coldata_degs$SampleID
coldata_degs$Status <- as.factor(coldata_degs$Status)
cellmeans_degs <- round(t(cellmeans_degs),0)
cellmeans_degs <- cellmeans_degs[, rownames(coldata_degs)]

dsd_degs <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(countData = cellmeans_degs, colData = coldata_degs, design = ~ Status))

dsd_degs <- suppressMessages(DESeq2::DESeq(dsd_degs, test = "LRT",reduced = ~ 1,fitType = 'parametric',
                                           sfType = 'poscounts',
                                           betaPrior = F))

res_degs <- as.data.frame(DESeq2::results(dsd_degs))
pvalues_degs <- as.numeric(res_degs$pvalue)

DE_genes <- res_degs %>% 
  mutate(avg_logFC = log(2^(log2FoldChange)), FC = 2^(log2FoldChange), FC = ifelse(FC<1, -1/FC, FC), gene = rownames(res_degs), 
         p_val = pvalue, p_val_adj = padj) %>%
  dplyr::select(p_val, FC, avg_logFC, p_val_adj, gene) %>%
  arrange(p_val)

write.csv(DE_genes, file.path(temp.dir, DEs_file),row.names=F)

