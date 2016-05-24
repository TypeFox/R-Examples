#'@keywords internal
#'@importFrom stats dgamma
logposterior_DPMST <- function(z, xi, psi, Sigma, df, B, hyper, c, m, alpha, n, a, b, diagVar){
  
  res <- NA
  
  indfull <- which(m!=0)
  mfull <- m[indfull]
  K <- length(indfull)
  
  if(!is.list(xi) && is.null(dim(xi))){
    log_vrais <- sum(mvstpdf(x = z, xi = xi, sigma = Sigma, psi = psi, df=df, Log=TRUE))
    if(!diagVar){
      log_prior_NNiW <-  sum(dNNiW(xi, psi, B, Sigma, hyperprior=hyper, log=TRUE))
    }else{
      betas <- diag(Sigma)
      beta0 <- diag(hyper$lambda)
      log_prior_NNiW <- sum(stats::dgamma(x=betas,
                                          shape=hyper$nu,
                                          rate=1/beta0, log=TRUE))
    }
  }else{
    
    if(!is.list(xi)){
      U_xi_full <- sapply(indfull, function(j) xi[, j])
      U_psi_full <- sapply(indfull, function(j) psi[, j])
      U_Sigma_full <- lapply(indfull, function(j) Sigma[, ,j])
      U_df_full <- sapply(indfull, function(j) df[j])
    }else{
      U_xi_full <- sapply(xi, "[")
      U_psi_full <- sapply(psi, "[")
      U_Sigma_full <- Sigma
      U_df_full <- sapply(df, "[")
    }
    if(nrow(z)==1){
      U_Sigma_full <- lapply(U_Sigma_full, FUN=matrix, nrow=1, ncol=1)
      U_xi_full <- matrix(U_xi_full, nrow=1)
      U_psi_full <- matrix(U_psi_full, nrow=1)
    }
    log_lik <- mvstlikC(x=z, c=c, clustval=indfull,
                        xi=U_xi_full, psi=U_psi_full, sigma=U_Sigma_full, df=U_df_full,
                        loglik=TRUE)
    log_vrais <- log_lik$total
    
    if(!diagVar){
      log_prior_NNiW <-  sum(dNNiW(lapply(indfull, function(j) xi[, j]),
                                   lapply(indfull, function(j) psi[, j]),
                                   lapply(indfull, function(j) B[, , j]),
                                   U_Sigma_full,
                                   hyperprior = hyper, log = TRUE))
    }else{
      betas <- lapply(U_Sigma_full, diag)
      beta0 <- diag(hyper$lambda)
      S <- lapply(betas, function(b){sum(stats::dgamma(x = b, shape = hyper$nu,
                                                       rate = 1/beta0, log = TRUE))})
      log_prior_NNiW <- sum(unlist(S))
    }
    
  }
  
  log_prior_alpha <- stats::dgamma(alpha, shape=a, scale=1/b, log=TRUE)
  
  log_clustering <- sum(c(lgamma(alpha), K*log(alpha), lgamma(mfull), -lgamma(alpha + n)))
  
  res <- c(log_vrais, log_clustering, log_prior_NNiW, log_prior_alpha)
  
  return(res)
  
}