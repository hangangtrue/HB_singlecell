rm(list = ls())
library(Seurat)
library(BiocParallel)
library(dplyr)


#-------------------------------------------------------------------------------
#                            main function
#-------------------------------------------------------------------------------

info_fcn <- function(indx, mat, prob, meta, pseudo_bulk, bulk_mean, bulk_se, niter = 100){
  ### generate result.tmp (22*7) in code.R
  gname = rownames(mat)[indx]
  sub <- meta[,3:4]
  sub <- sub[!duplicated(sub[,1]),]
  subject <- unique(sub$Donor)
  
  rownames(pseudo_bulk)<-pseudo_bulk$X
  pseudo_bulk<-pseudo_bulk[,-1]
  
  ind_data <- matrix(NA,nrow = length(subject), ncol = 5)
  colnames(ind_data)<-c("Pheno","Num_Cells","Prob_mean","Negative_log","Mean")
  for(j in 1:length(subject)){
    
    cell<-which(meta$Donor==subject[j])
    s<-which(sub$Donor==subject[j])
    
    ind_data[j,"Pheno"]<- as.numeric(sub$pheno[s]=="IPF")    
    ind_data[j,"Num_Cells"]<- length(which(meta$Donor==subject[j]))
    ind_data[j,"Prob_mean"]<- prob[which(prob$Donor==subject[j]),2]
    ind_data[j,"Negative_log"]<- -log(ind_data[j,3])
    ind_data[j,"Mean"]<- pseudo_bulk[indx, which(colnames(pseudo_bulk)==subject[j])]
  }
  
  rownames(ind_data)<-subject
  
  # design matrix from ind_data
  Y <- ind_data[,"Mean"]
  X <- matrix(c(rep(1, length(Y)), ind_data[,"Pheno"], ind_data[,"Negative_log"]),ncol=3)
  
  
  # frequentist estimates
  freq_mu1 <- solve(t(X)%*%X)%*%(t(X)%*%Y)
  sigmasq <- c(var(Y-X%*%freq_mu1))  # estimated sigma
  freq_sig1 <- solve(t(X) %*% X) * sigmasq
  
  # Prior information
  nonin_mu1.string <- c(0,0,0)
  nonin_mu1 <- matrix(nonin_mu1.string,nrow=3)
  nonin_sig1 <- matrix(0,3,3)
  diag(nonin_sig1) <- 100
  
  ### Input informative prior
  nonin_mu1[2,1] <- bulk_mean[indx]
  nonin_sig1[2,2] <- (bulk_se[indx])^2

  
  # Posterior mean and variance
  nonin_mu1_pos <- solve(solve(nonin_sig1)+t(X)%*%X)%*%(solve(nonin_sig1)%*%nonin_mu1+t(X)%*%Y)
  nonin_sig1_pos <- solve(solve(nonin_sig1)+t(X)%*%X)*sigmasq
  
  # Bayesian frequentist hybrid inference
  # The only Bayesian parameter is pheno
  hb_mu <-  bulk_mean[indx]
  hb_sig <- (bulk_se[indx])^2

  # Step 1. Initialize parameters
  # Get initial value from frequentist method
  freq_mu1 <- solve(t(X)%*%X)%*%(t(X)%*%Y)
  # Step 2-3. At the current value of Beta, draw the posterior and estimate the
  # Bayesian parameter: beta1
  # frequentist parameter: beta0, beta2
  cols <- c(1,3)
  Y_temp <- Y-X[,cols]%*%freq_mu1[cols,]
  X_temp <- X[,2]
  hb_mu1_pos <- (solve(solve(hb_sig)+t(X_temp) %*% X_temp)) %*% (solve(hb_sig) %*% hb_mu+t(X_temp) %*% Y_temp)
  # Step 4. Given the estimate in Step 3, estimate/update the frequentist parameter
  Y_temp <- Y-(X[,2])*as.vector(hb_mu1_pos)
  X_temp <- X[,cols]
  hb_mu1_fre <- solve(t(X_temp) %*% X_temp) %*% (t(X_temp) %*% Y_temp)
  # Step 5. Iterate steps 2-4 (for example, 100 times, as in EM algorithm)
  for (i in seq_len(niter)){
    Y_temp <- Y-X[,cols]%*%hb_mu1_fre
    X_temp <- X[,2]
    hb_mu1_pos <- solve(solve(hb_sig)+t(X_temp)%*%X_temp)%*%(solve(hb_sig)*hb_mu+t(X_temp)%*%Y_temp)
    Y_temp <- Y-X[,2]%*%hb_mu1_pos
    X_temp <- X[,cols]
    hb_mu1_fre <- solve(t(X_temp)%*%X_temp)%*%(t(X_temp)%*%Y_temp)
  }
  hb_est <- c(hb_mu1_fre[1,1],hb_mu1_pos,hb_mu1_fre[2,1])  ## Stored as a vector, not a matrix
  
  EST <- c(freq_mu1, nonin_mu1_pos, hb_est)
  EST <- matrix(EST,nrow=3)
  
  # Variance estimate
  # Frequentist
  var_freq <- diag(freq_sig1)
  
  # Bayesian
  var_bayes <- diag(nonin_sig1_pos)
  
  # Hybrid Bayesian
  sigmasq_hbf <- var(Y_temp-X_temp%*%hb_mu1_fre)
  hb_freq_sig <- solve(t(X_temp)%*%X_temp)*as.vector(sigmasq_hbf)
  
  Y_temp <- Y-X[,cols]%*%hb_mu1_fre
  X_temp <- X[,2]
  sigmasq_hbb <- var(Y_temp-X_temp%*%hb_mu1_pos)
  hb_bayes_sig_pos <- solve(solve(hb_sig)+t(X_temp)%*%X_temp)%*%sigmasq_hbb
  
  var_hb <- c(hb_freq_sig[1,1], hb_bayes_sig_pos, hb_freq_sig[2,2])
  
  VAR <- c(var_freq, var_bayes, var_hb)
  VAR <- matrix(VAR, nrow=3)
  
  STD <- sqrt(VAR)
  
  # Using est and std to calculate p-value
  z_score <- matrix(0,nrow=nrow(EST),ncol=ncol(EST))
  p_value <- matrix(0,nrow=nrow(EST),ncol=ncol(EST))
  
  for (i in 1:nrow(EST)){
    for (j in 1:ncol(EST)){
      z_score[i,j] <- EST[i,j]/STD[i,j]
      if (pnorm(z_score[i,j])<0.5) 
      {p_value[i,j] <- 2*pnorm(z_score[i,j])}
      else   
      {p_value[i,j] <- 2*(1-pnorm(z_score[i,j]))}
    }
  }
  saveRDS(list(EST = EST, STD = STD, pval = p_value), paste0("d230410_info_output_", gname, ".rds"))
  
}

#-------------------------------------------------------------------------------
#                            load data
#-------------------------------------------------------------------------------


load("CountMatrix.Rdata") # name: matrix
# matrix = matrix[1:10,]
indx_vec = seq_len(nrow(matrix)) # number of genes

bulk = read.csv("Bulk_IPF_GSE150910.csv")
bulk <- bulk %>% 
  mutate(bulk_indx = seq_len(nrow(bulk)))

matrix_gname_df <- data.frame(geneName = rownames(matrix), indx = seq_len(nrow(matrix)))
bulk_match_matrix <- matrix_gname_df %>% left_join(bulk, by = "geneName")
saveRDS(bulk_match_matrix, "bulk_match_matrix.rds")

#-------------------------------------------------------------------------------
#                            load data
#-------------------------------------------------------------------------------
Start = Sys.time()
BiocParallel::bplapply(
  X = indx_vec, 
  FUN = info_fcn, 
  BPPARAM = MulticoreParam(5),
  mat = matrix, 
  prob = read.csv("Prob_mean.csv",header=TRUE),
  meta = read.csv("GSE135893_meta.csv", header=TRUE),
  pseudo_bulk = read.csv("ind_data_1.csv",header=TRUE), 
  bulk_mean = bulk_match_matrix$diff,
  bulk_se = bulk_match_matrix$SE)
End = Sys.time()
End - Start
# Time difference of 3.916346 mins



