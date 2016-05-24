clime <- function(x, lambda=NULL,
                  nlambda=ifelse(is.null(lambda),100,length(lambda)),
                  lambda.max=0.8,
                  lambda.min=ifelse(nrow(x)>ncol(x), 1e-4, 1e-2),
                  sigma=FALSE,
                  perturb=TRUE,
                  standardize=TRUE,
                  logspaced=TRUE,
                  linsolver=c("primaldual", "simplex"),
                  pdtol=1e-3, pdmaxiter=50
                  )
{
  lpfun <- match.arg(linsolver, c("primaldual", "simplex"))
  
  if (sigma) {
    if (is.matrix(x)) {
      Sigma <- x
    } else {
      Sigma <-as.matrix(x)
    }
    p <- ncol(Sigma)
    x <- NULL
  } else {
    n <- nrow(x)
    p <- ncol(x)
    
    if (is.null(lambda)) {
      if (logspaced) {
        lambda <- 10^(seq(log10(lambda.min), log10(lambda.max), length.out=nlambda))
      } else {
        lambda <- seq(lambda.min, lambda.max, length.out=nlambda)
      }
    }
    

    if (standardize)  x <- scale(x)
    Sigma <- cov(x)*(1-1/n)
  }
  
  ## Set to perturbed Sigma to have conditional number p
  eigvals <- eigen(Sigma, only.values=T)$values
  if (is.logical(perturb)) {
      if (perturb) { 
          perturb <- max(max(eigvals) - p*min(eigvals), 0)/(p-1)
      } else {
          perturb <- 0
      }
  }
  
  Sigma <- Sigma+diag(p)*perturb
  emat <- diag(p)

  Omegalist <- vector("list", nlambda)
  if (lpfun == "simplex") {
    for (jl in 1:nlambda) {
      Omega <- matrix(0, nrow=p, ncol=p)
      lam <- lambda[jl]
      for (j in 1:p) {
        beta <- linprogS(Sigma, emat[,j], lam)
        Omega[,j] <- beta
      }
      Omegalist[[jl]] <- Omega*(abs(Omega)<= abs(t(Omega)))+ t(Omega)*(abs(Omega)> abs(t(Omega)))
    }
  }

  if (lpfun == "primaldual") {
    Omega0 <- solve(Sigma)
    
    for (jl in 1:nlambda) {
      Omega <- matrix(0, nrow=p, ncol=p)
      lam <- lambda[jl]
      for (j in 1:p) {
        beta <- linprogPD(Omega0[,j], Sigma, emat[,j], lam, pdtol, pdmaxiter)
        Omega[,j] <- beta
      }
      Omegalist[[jl]] <- Omega*(abs(Omega)<= abs(t(Omega)))+ t(Omega)*(abs(Omega)> abs(t(Omega)))
    }
  }

  
  

  outlist <- list(Omegalist=Omegalist, x = x, lambda = lambda, perturb=perturb, standardize = standardize, lpfun=lpfun)
  class(outlist) <- c("clime")
  return(outlist)
}

