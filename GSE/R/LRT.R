###################################################################
## This file contains functions extracted from robmiss packages
## written by Mike Danilov
slrt <- function(S, trueS){
    ## Standardized LRT-statistic for testing if cov S is the same as trueS
    if (missing(trueS)) trueS <- diag(rep(1,ncol(S)))
    if (is.character(trueS)) trueS <- .std.cov(p=ncol(S),trueS)
    if (!is.matrix(trueS)) stop("True covariance (trueS) misspecified")

    return(.lrt(solve(trueS)%*%S))
}
  
.lrt <- function(S)
{
  ## LRT-statistic for testing if covariance is identity  
  if (sum(is.na(S))>0) return(NA)
  return(sum(diag(S))-log(det(S))-ncol(S))
}
  
.std.cov <- function(p, ct, seed=2727)
{
  if (!exists(".Random.seed")) runif(1) # to initialize .Random.seed
  seed.before <- .Random.seed
  if (!is.null(seed)) set.seed(seed)
  res <- switch(ct,
                h=cor(mvrnorm(p+1, mu=rep(0,p), Sigma=diag(rep(1,p)))),
                m=cor(mvrnorm(round(p+5), mu=rep(0,p), Sigma=diag(rep(1,p)))),
                l=cor(mvrnorm(p*10, mu=rep(0,p), Sigma=diag(rep(1,p)))),
                i=diag(rep(1,p)))
  .Random.seed <- seed.before
  return(res)
}

