computePCC <- function(x)
{
  x <- -cov2cor(x)
  diag(x) <- 0
  x <- as.matrix(forceSymmetric(x))
  return(x)
}

computePDC <- function(beta,kappa){
  sigma <- solve(kappa)
  t(beta / sqrt(diag(sigma) %o% diag(kappa) + beta^2))
}

graphicalVAR <-
function(
  data, # A n by p data frame containing repeated measures
  nLambda = 50, # Either single value or vector of two corresponding to c(kappa, beta)
  verbose = TRUE,
  gamma = 0.5,
  lambda_beta,
  lambda_kappa, maxit.in = 100, maxit.out = 100
  ){
  
  # Check input:
  if (is.data.frame(data)){
    data <- as.matrix(data)
  }
 
  stopifnot(is.matrix(data))
  
  

  
  Nvar <- ncol(data)
  Ntime <- nrow(data)

  # Center data:
  data <- scale(data, TRUE, FALSE)
  
  # Compute current and lagged data:
  data_c <- data[-1,,drop=FALSE]
  data_l <- data[-nrow(data),,drop=FALSE]
  
  # Generate lambdas (from SparseTSCGM package):
  if (missing(lambda_beta) | missing(lambda_kappa)){
    lams <- SparseTSCGM_lambdas(data_l, data_c, nLambda)
    if (missing(lambda_beta)){
      lambda_beta <- lams$lambda_beta
    }
    if (missing(lambda_kappa)){
      lambda_kappa <- lams$lambda_kappa
    }
  }
  
  Nlambda_beta <- length(lambda_beta)
  Nlambda_kappa <- length(lambda_kappa)
  
  
  # Expand lambda grid:
  lambdas <- expand.grid(kappa = lambda_kappa, beta = lambda_beta)
  Estimates <- vector("list", nrow(lambdas))
  
  ### Algorithm 2 of Rothmana, Levinaa & Ji Zhua
  if (verbose){
    pb <- txtProgressBar(0, nrow(lambdas), style = 3) 
  }
  for (i in seq_len(nrow(lambdas))){
    Estimates[[i]] <- Rothmana(data_l, data_c, lambdas$beta[i],lambdas$kappa[i], gamma=gamma,maxit.in=maxit.in, maxit.out = maxit.out)
   if (verbose){
     setTxtProgressBar(pb, i)
   } 
  }
  if (verbose){
    close(pb)
  }
#   
#   logandbic <- LogLik_and_BIC(data_l, data_c, Estimates)
#   lambdas$bic <- logandbic$BIC
#   lambdas$loglik <- logandbic$logLik
  lambdas$ebic <- sapply(Estimates,'[[','EBIC')
  # Which minimal BIC:
  min <- which.min(lambdas$ebic)
  Results <- Estimates[[min]]

  # Standardize matrices (Wild et al. 2010)
  # partial contemporaneous correlation (PCC) 
  Results$PCC <- computePCC(Results$kappa)
  Results$PDC <- computePDC(Results$beta, Results$kappa)  

  Results$path <- lambdas
  Results$labels <- colnames(data)

  colnames(Results$beta) <- rownames(Results$beta) <- colnames(Results$kappa) <- rownames(Results$kappa) <-
  colnames(Results$PCC) <- rownames(Results$PCC) <- colnames(Results$PDC) <- rownames(Results$PDC) <-
  Results$labels
Results$gamma <- gamma

  class(Results) <- "graphicalVAR"
  
  return(Results)
}
