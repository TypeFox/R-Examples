"float" <-
function(object, factor, iter.max = 50)
{
  float.variance <- function(V, tol=1.0e-3, iter.max = 50)
    {
      ## Calculate floated variances for variance matrix V, which is
      ## assumed to represent a set of treatment contrasts
    
      m <- nrow(V)
      if (!is.matrix(V) || ncol(V) != m || m == 1)
        stop ("V must be a square matrix of size 2 x 2 or more")

      evals <- eigen(V, only.values=TRUE)$values
      if(any(evals < 0))
        stop("V not positive definite")
    
      ## Starting values from Easton et al (1991)
      R <- V - diag(diag(V))
      V00 <- sum(R)/(m * (m-1))
      V10 <- apply(R, 1, sum)/(m-1)
      fv <- c(V00, V00 - 2 * V10 + diag(V))
    
      for(iter in 1:iter.max) {
        w <- 1/fv
        S <- sum(w)
        w1 <- w[-1]/S
        ##Augment data matrix
        V10 <- as.vector(V %*% w1)
        V00 <- as.vector(1/S + t(w1) %*% V %*% w1)
        ##Calculate new estimates
        fv.old <- fv
        fv <- c(V00, V00 - 2 * V10 + diag(V))
        ## Check convergence
        if(max(abs(fv.old - fv)/fv) < tol)
          break
      }
      if (iter == iter.max)
        warning("Floated variance estimates did not converge")
  
      Vmodel.inv <- S * (diag(w1) - w1 %*% t(w1))
      evals <- 1/(eigen(V %*% Vmodel.inv, only.values=TRUE)$values)
      divergence <- sum(1/evals - 1 + log(evals))/2
      return(list(variance=fv, error.limits=sqrt(range(evals)),
                  divergence=divergence))
    }

  if (is.null(object$xlevels)) {
    stop("No factors in model")
  }
  
  if (missing(factor)) {
    i <- 1
    factor <- names(object$xlevels)[1]
  }
  else {
    i <- pmatch(factor, names(object$xlevels))
    if (is.na(i)) {
      stop(paste("Factor",i,"not found in model"))
    }
  }
    
  xcontrasts <- object$contrasts[[i]]
  xlevels <- object$xlevels[[i]]
  xname <- names(object$xlevels)[i]
  nlevels <- length(xlevels)

  ## Extract the coefficients and variance matrix for a single factor
  ## from object
        
  if (nlevels <= 2) {
    stop ("Floated variances undefined for factors with less than 3 levels")
  }
        
  ## Get contrast matrix
  C <- if (is.matrix(xcontrasts)) {
    xcontrasts
  }
  else {
    get(xcontrasts, mode="function")(xlevels)
  }
  if (qr(C)$rank < nlevels - 1) {
    stop ("Impossible to reconstruct treatment contrasts")
  }

  ## Get coefficients and variance matrix
  if(is.null(cnames <- colnames(C)))
    cnames <- 1:(nlevels-1)
  contr.names <- paste(xname, cnames, sep="")
  coef <- coef(object)[contr.names]
  V <- vcov(object)[contr.names, contr.names]
        
  ## Convert to treatment contrast parameterization
  if (identical(xcontrasts, "contr.treatment")) {
    V.tc <- V
    coef.tc <- c(0, coef)
  }
  else {
    D.inv <- cbind(rep(-1,nlevels-1), diag(nlevels-1))
    S <- D.inv %*% cbind(rep(1, nlevels), C) 
    S <- S[,-1]
    ## coefficients
    coef.tc <- c(0, S %*% coef)
    ## If we find a baseline level (implicitly defined
    ## by having a row of zeros in the contrast matrix)
    ## then adjust the coefficients
    is.base <- apply(abs(C), 1, sum) == 0
    if (any(is.base))
      coef.tc <- coef.tc - coef.tc[is.base]
    ## variance matrix
    V.tc <- S %*% V %*% t(S)
  }
  names(coef.tc) <- xlevels

  float.out <- float.variance(V.tc, iter.max = iter.max)
  var <- float.out$var
  names(var) <- xlevels
  ans <- list(coef=coef.tc, var=var, limits=float.out$error.limits,
              factor=factor)
  class(ans) <- "floated"
  return(ans)
}
