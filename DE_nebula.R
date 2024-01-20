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
write.csv(s_g, file = file.path("nebula",selected_tDE), row.names = FALSE)

### Add effect size (fold change) to the selected genes in the Disease group
combined[["RNA"]]@counts[rownames(combined) %in% s_g, combined$status_p == "Disease"] <- combined[["RNA"]]@counts[rownames(combined) %in% s_g, combined$status_p == "Disease"]*f

### Assign groups
combined$label <- combined$status_p
DefaultAssay(combined) <- "RNA"
meta_group = levels(combined$status_p)
ident.1 = meta_group[1]
ident.2 = meta_group[2]
Idents(combined) <- combined$C2_named

#### Nebula
pct <- 0
min.cells.feature = 0
min.cells.group = 0
logfc.threshold = 0
test.use <- "nebula"
i <- "C2-1: Neurons"
DEs_file <- paste("Mice_", m, "_Cells_", ce, "_FC_", f, "_Rep_",j,"_data.csv", sep = "")
DefaultAssay(combined) <- "RNA"

### DE genes
sub_grp <- combined %>%
  subset( (Idents(combined) == i) & (status_p %in% c(ident.1, ident.2)))
counts <- sub_grp[["RNA"]]@counts
orig.ident <- as.character(sub_grp$Sample_ID)
grp <- data.frame(group = sub_grp$status_p)
df = model.matrix(~group, data=grp)
results <- nebula(counts, orig.ident, pred=df)

pct.1 <- rowMeans(counts[, grp$group == ident.1] > 0)
pct.2 <- rowMeans(counts[, grp$group == ident.2] > 0)
pct_df <- data.frame(pct.1 = pct.1, pct.2 = pct.2, gene = names(pct.1))

DE_genes <- results$summary %>%
  mutate(FC = exp(results$summary[,2]), FC = ifelse(FC<1, -1/FC, FC),
         p_val = results$summary[,6], p_val_adj = p.adjust(p_val, method = "fdr"),
         cluster = i, overdispersion_sub = results$overdispersion[,1],
         overdispersion_cell = results$overdispersion[,2], convergence = results$convergence,
         convergence = ifelse(convergence == 1, "Yes", "No"))

DE_genes <- merge(DE_genes, pct_df, by = "gene") %>%
  dplyr::select(p_val, FC, pct.1, pct.2, overdispersion_sub, overdispersion_cell, convergence, p_val_adj, cluster,
                gene) %>%
  arrange(p_val)

write.csv(DE_genes, file.path(temp.dir, DEs_file),row.names=F)
