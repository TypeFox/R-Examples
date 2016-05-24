ridge.proj <- function(x, y,
                       family = "gaussian",
                       standardize = TRUE,
                       lambda = 1,
                       betainit = "cv lasso",
                       sigma = NULL,
                       suppress.grouptesting = FALSE,
                       multiplecorr.method = "holm", N = 10000)
{
  ## Purpose:
  ## calculation of the ridge projection method proposed in
  ## http://arxiv.org/abs/1202.1377 P.Buehlmann
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## x: the design matrix
  ## y: the response vector
  ## N: number of simulations for the WY procedure
  ## ----------------------------------------------------------------------
  ## Return values:
  ## pval: the individual testing p-values for each parameter
  ## pval.corr:  the multiple testing corrected p-values for each parameter
  ## betahat:    the initial estimate by the scaled lasso of \beta^0
  ## bhat:       the de-sparsified \beta^0 estimate used for p-value calculation
  ## sigmahat:   the \sigma estimate coming from the scaled lasso
  ## ----------------------------------------------------------------------
  ## Author: Peter Buehlmann (initial version),
  ##         adaptations by L. Meier and R. Dezeure

  n <- nrow(x)
  p <- ncol(x)

  if(standardize)
    sds <- apply(x, 2, sd)
  else
    sds <- rep(1, p)

  ## *center* (scale) the columns
  x <- scale(x, center = TRUE, scale = standardize)

  dataset <- switch(family,
                    "gaussian" = {
                      list(x = x, y = y)
                    },
                    {
                      switch.family(x = x, y = y,
                                    family = family)
                    })

  x <- scale(dataset$x, center = TRUE, scale = FALSE)
  y <- scale(dataset$y, scale = FALSE)
  y <- as.numeric(y)

  ## Warning: should we allow user to still specify their
  ## own Z here?
  ## Z <- NULL##will have to recalculate Z
  
  ## force sigmahat to 1 when doing glm!
  if(family == "binomial")
    sigma <- 1

  ## these are some old arguments that are still used in the code below
  ridge.unprojected <- TRUE

  ## *center* (scale) the columns
  x <- scale(x, center = TRUE, scale = standardize) 
  y <- scale(y, scale = FALSE)
  
  y  <- as.numeric(y)

  biascorr <- Delta <- numeric(p)
  
  ##lambda <- 1  ## other choices?

  h1 <- svd(x)
  
  ## Determine rank of design matrix
  ## Overwrite h1 with the version corresponding to non-zero singular values
  sval  <- h1$d
  tol   <- min(n, p) * .Machine$double.eps
  rankX <- sum(sval >= tol * sval[1])

  h1$u <- h1$u[,1:rankX]
  h1$d <- h1$d[1:rankX]
  h1$v <- h1$v[,1:rankX]
  ## End overwrite h1 with the version corresponding to non-zero singular values

  Px               <- tcrossprod(h1$v)
  Px.offdiag       <- Px
  diag(Px.offdiag) <- 0 ## set all diagonal elements to zero

  ## Use svd for getting the inverse for the ridge problem
  hh <- h1$v %*% ((h1$d / (h1$d^2 + lambda)) * t(h1$u))
  
  ## Ruben Note: here the h1$d^2/n 1/n factor has moved out, the value of lambda
  ## used is 1/n! See also comparing to my version
  
  cov2      <- tcrossprod(hh) 
  diag.cov2 <- diag(cov2)

  ## Estimate regression parameters (initial estimator) and noise level sigma
  initial.estimate <- initial.estimator(betainit = betainit, sigma = sigma,
                                        x = x,y = y)
  beta.lasso <- initial.estimate$beta.lasso
  hat.sigma2 <- initial.estimate$sigmahat^2
  
  ## bias correction
  biascorr <- crossprod(Px.offdiag, beta.lasso)

  ## ridge estimator
  hat.beta <- hh %*% y

  hat.betacorr <- hat.beta - biascorr

  if(ridge.unprojected){ ## bring it back to the original scale
    hat.betacorr <- hat.betacorr / diag(Px)
  }
  
  ## Ruben Note: a_n = 1 / scale.vec, there is no factor sqrt(n) because this
  ## falls away with the way diag.cov2 is calculated see paper
  scale.vec  <- sqrt(hat.sigma2 * diag.cov2)
  
  ## 1, is coupled with 2!! don't shut this off and not the other or vice versa
  if(ridge.unprojected){
    scale.vec <- scale.vec / abs(diag(Px))
  }
  
  hat.betast <- hat.betacorr / scale.vec ## sqrt(hat.sigma2 * diag(cov2))
  Delta      <- apply(abs(Px.offdiag), 2, max) * (log(p) / n)^(0.45) / scale.vec
  ## new
  if(ridge.unprojected ){ ## 2
    Delta <- Delta / abs(diag(Px))
    ## to put it on the same scale when scale.vec is different!
  }
     
  hat.gamma1 <- abs(hat.betast)

  #########################
  ## Individual p-values ##
  #########################
  
  res.pval <- c(pmin(2 * pnorm(hat.gamma1 - Delta, lower.tail = FALSE), 1))

  #########################################
  ## Multiple testing corrected p-values ##
  #########################################
  
  if(multiplecorr.method == "WY"){
    ## Westfall-young like procedure 
    pcorr <- p.adjust.wy(cov = cov2, pval = res.pval, N = N)
  }else{
    if(multiplecorr.method %in% p.adjust.methods){
      pcorr <- p.adjust(res.pval, method = multiplecorr.method)
    }else{
      stop("Unknown multiple correction method specified")
    }
  }

  ##########################
  ## Confidence Intervals ##
  ##########################

  scaleb <- 1 / scale.vec
  se     <- 1 / scaleb

  ## need to multiply Delta with se because it is on the scale of
  ## standard normal dist and we want to bring it to the distribution of
  ## \hat{\beta}_j

  ##############################################
  ## Function to calculate p-value for groups ##
  ##############################################
  if(suppress.grouptesting){
    group.testing.function <- NULL
    cluster.group.testing.function <- NULL
  }else{
    pre <- preprocess.group.testing(N = N, cov = cov2, conservative = FALSE)
    
    group.testing.function <- function(group, conservative = FALSE){
      calculate.pvalue.for.group(brescaled    = hat.betast,
                                 group        = group,
                                 individual   = res.pval,
                                 Delta        = Delta,
                                 ##correct    = TRUE,
                                 conservative = conservative,
                                 zz2          = pre)
    }
    
    cluster.group.testing.function <-
      get.clusterGroupTest.function(group.testing.function=group.testing.function,
                                    x = x)
  }
  
  out <- list(pval          = res.pval,
              pval.corr     = pcorr,
              groupTest     = group.testing.function,
              clusterGroupTest = cluster.group.testing.function,
              sigmahat      = sqrt(hat.sigma2),
              standardize   = standardize,
              sds           = sds,
              bhatuncorr    = hat.beta / sds,
              biascorr      = biascorr / sds,
              normalisation = 1 / scale.vec,
              bhat          = hat.betacorr / sds,
              se            = se / sds,
              delta         = Delta,
              betahat       = beta.lasso / sds,
              family        = family,
              method        = "ridge.proj",
              call          = match.call())
  
  names(out$pval) <- names(out$pval.corr) <- names(out$bhat) <-
    names(out$sds) <- names(out$se) <- names(out$delta) <- names(out$betahat) <-
      colnames(x)

  class(out) <- "hdi"

  return(out)
}
