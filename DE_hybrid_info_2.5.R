library(hierarchicell)
library(nebula)
library(Libra)
library(pbmcapply)

source("main_fn.R")


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


######################################################################
############      Manipulate data to generate setup      #############
######################################################################
Hypo <- "Neuro_normal_chow.RData"
load(file.path(Hypo))
#### The genes are always true DE genes, just to match the previous practice. It's OK to include them every iteration. 
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
write.csv(s_g, file = file.path("hybrid_info_2.5",selected_tDE), row.names = FALSE)

### Add effect size (fold change) to the selected genes in the Disease group
combined[["RNA"]]@counts[rownames(combined) %in% s_g, combined$status_p == "Disease"] <- combined[["RNA"]]@counts[rownames(combined) %in% s_g, combined$status_p == "Disease"]*f

### Assign groups
combined$label <- combined$status_p
DefaultAssay(combined) <- "RNA"
meta_group = levels(combined$status_p)
ident.1 = meta_group[1]
ident.2 = meta_group[2]
Idents(combined) <- combined$C2_named
cores = 5

### Calculate and store FC
counts_mt <- as.matrix(combined[["RNA"]]@counts)
mu.1 <- rowMeans(counts_mt[, combined$label == ident.1])
mu.2 <- rowMeans(counts_mt[, combined$label == ident.2])
FC <- data.frame(gene=rownames(counts_mt), FC = mu.1/mu.2)

data_file <- paste("Mice_", m, "_Cells_", ce, "_FC_", f, "_Rep_",j,"_afc.csv", sep = "")
write.csv(FC, file.path(temp.dir, data_file),row.names=F)


######################################################################
##############                hybrid_info             ################
######################################################################
###### Prepare pseudo bulk and meta data for hybrid
i <- "C2-1: Neurons"
test.use <- "hybrid_info"
DEs_file <- paste("Mice_", m, "_Cells_", ce, "_FC_", f, "_Rep_",j,"_data.csv", sep = "")

combined$Sample_ID <- gsub("_", "-", combined$Sample_ID)
meta <- data.frame(DonorID = paste(combined$status_p, combined$Sample_ID, sep = "_"))
genecounts <- t(as.matrix(combined@assays$RNA@counts))
genecounts <- cbind(meta,genecounts)
computemeans <- function(x){tapply(x,genecounts[,1],mean)}
cellmeans <- sapply(genecounts[,-1],computemeans)
## Input is a gene by subject matrix
pseudo_bulk <- t(cellmeans)
meta <- sapply(strsplit(colnames(pseudo_bulk), split = "_"), function(inner_list) inner_list[[1]])
meta[meta=="Disease"] <- 1
meta[meta=="HC"] <- 0
meta <- as.numeric(meta)
indx_vec <- seq_len(nrow(pseudo_bulk))

## Info prior for true DEs, noninfo for Null genes
Start = Sys.time()
DE_genes <- c()
for (g in 1:nrow(pseudo_bulk)) {
  gname = rownames(pseudo_bulk)[g]
  if(gname %in% s_g){
    results <- info_fcn(indx = g, meta, pseudo_bulk, prior_mean=2.5, prior_se=1, niter = 100)
    DE_genes <- rbind(DE_genes, results)
  } else {
    results <- info_fcn(indx = g, meta, pseudo_bulk, prior_mean=0, prior_se=1000, niter = 100)
    DE_genes <- rbind(DE_genes, results)
  }
  print(g)
}
End = Sys.time()
End - Start
# Time difference of 4.798463 mins
DE_genes <- DE_genes %>% 
  arrange(p_val)

## Writing results
write.csv(DE_genes, file.path(temp.dir, DEs_file),row.names=F)
