library(Seurat)
library(BiocParallel)
library(dplyr)


#-------------------------------------------------------------------------------
#                            main function
#-------------------------------------------------------------------------------

info_fcn <- function(indx, meta, pseudo_bulk, prior_mean, prior_se, niter = 100){

  gname = rownames(pseudo_bulk)[indx]
  
  ind_data <- matrix(NA,nrow = length(meta), ncol = 2)
  colnames(ind_data)<-c("Pheno","Mean")
  ind_data[,"Pheno"]<- meta   
  ind_data[,"Mean"]<- pseudo_bulk[indx, ]
  
  rownames(ind_data)<- colnames(pseudo_bulk)
  
  # design matrix from ind_data
  Y <- ind_data[,"Mean"]
  X <- matrix(c(rep(1, length(Y)), ind_data[,"Pheno"]),ncol=2)

  # frequentist estimates
  freq_mu1 <- solve(t(X)%*%X)%*%(t(X)%*%Y)
  sigmasq <- c(var(Y-X%*%freq_mu1))  # estimated sigma
  freq_sig1 <- solve(t(X) %*% X) * sigmasq
  
  # Prior information
  nonin_mu1.string <- c(0,0)
  nonin_mu1 <- matrix(nonin_mu1.string,nrow=2)
  nonin_sig1 <- matrix(0,2,2)
  diag(nonin_sig1) <- 100
  
  ### Input informative prior
  nonin_mu1[2,1] <- prior_mean
  nonin_sig1[2,2] <- (prior_se)^2
  # nonin_sig1[2,2] <- bulk_se[indx]
  
  
  # Posterior mean and variance
  nonin_mu1_pos <- solve(solve(nonin_sig1)+t(X)%*%X)%*%(solve(nonin_sig1)%*%nonin_mu1+t(X)%*%Y)
  nonin_sig1_pos <- solve(solve(nonin_sig1)+t(X)%*%X)*sigmasq
  
  # Bayesian frequentist hybrid inference
  # The only Bayesian parameter is pheno
  hb_mu <-  prior_mean
  hb_sig <- (prior_se)^2
  # hb_sig <- bulk_se[indx]
  
  # Step 1. Initialize parameters
  # Get initial value from frequentist method
  freq_mu1 <- solve(t(X)%*%X)%*%(t(X)%*%Y)
  # Step 2-3. At the current value of Beta, draw the posterior and estimate the
  # Bayesian parameter: beta1
  # frequentist parameter: beta0
  cols <- c(1)
  Y_temp <- Y-X[,cols]*freq_mu1[1]
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
  hb_est <- c(hb_mu1_fre[1,1],hb_mu1_pos)  ## Stored as a vector, not a matrix
  
  EST <- c(hb_est)
  EST <- matrix(EST)
  
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
  
  var_hb <- c(hb_freq_sig[1,1], hb_bayes_sig_pos)
  
  VAR <- c(var_hb)
  VAR <- matrix(VAR)
  
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
  results <- data.frame(gene = gname, EST = EST, STD = STD, p_val = p_value)
  results <- results[2,]
  
  return(results)
 
}
