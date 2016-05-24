Mort1Dsmooth_optimize <-
function(x, y, offset, wei, psi2,
                                  B, DtD,
                                  a.init,
                                  MON, TOL1, TOL2,
                                  RANGE, MAX.IT, MET){
  ## the function first applies a rough cleversearch
  ## then around the optimal values,
  ## it search with a fine cleversearch
  ## a more precise value for the smoothing parameter

  ## Input:
  ## x: abcissae of data
  ## y: count response
  ## offset: an a priori known component
  ## wei: weigths
  ## psi2: overdispersion parameter
  ## B: B-splines basis
  ## DtD: inner product of the matrix of differences
  ## a.init: starting coefficients
  ## MON: logical on monitoring
  ## TOL1: relative convergence tolerance
  ## TOL2: difference between two adjacent
  ##       smoothing parameters in the second grid search,
  ##       log-scale.
  ## RANGE: range in between grid-search
  ## MAX.IT: the maximum number of iterations
  ## MET: criterion for optimizing lambda (bic or aic)
  
  ## Output: a list containing
  ## lambda: the optimal smoothing parameter

  ## object function 
  Mort1Dsmooth_opt_ic <- function(X){
    FIT <- Mort1Dsmooth_estimate(x = x, y = y,
                                 offset = offset, 
                                 wei = wei,
                                 psi2 = psi2, B = B,
                                 lambda = X, DtD = DtD, 
                                 a.init = a.init,
                                 MON=MON, TOL1=TOL1,
                                 MAX.IT=MAX.IT)
    return(ifelse(MET==2, FIT$aic, FIT$bic))
  }

  
  ## ROUGH SEARCHING
  ## possible lambdas
  lambdas.0 <- seq(log10(RANGE[1]),log10(RANGE[2]),
                   TOL2*4)
  ## length of possible lambdas
  by.lambda.0 <- length(lambdas.0)
  ## starting values for lambdas in cleversearch
  l.st.0 <- median(lambdas.0)
  ## optimizing lambda
  lambda.hat0 <- cleversearch(fn=Mort1Dsmooth_opt_ic,
                              lower=log10(RANGE[1]),
                              upper=log10(RANGE[2]),
                              startvalue=l.st.0,
                              ngrid=by.lambda.0,
                              logscale=TRUE,
                              verbose=FALSE)[[1]]
  ## FINE SEARCHING
  ## sub-optimal lambda=starting value in fine cleversearch
  l.st <- log10(lambda.hat0)
  ## range in which search
  min.l <- max(l.st-TOL2*4, log10(RANGE[1]))
  max.l <- min(l.st+TOL2*4, log10(RANGE[2]))
  ## possible lambdas
  lambdas <- seq(min.l, max.l, TOL2)
  ## length of possible lambdas
  by.lambda <- length(lambdas)
  ## optimizing lambda
  lambda.hat <- cleversearch(fn=Mort1Dsmooth_opt_ic,
                             lower=min.l,
                             upper=max.l,
                             startvalue=l.st,
                             ngrid=by.lambda,
                             logscale=TRUE,
                             verbose=FALSE)[[1]]

  return(lambda.hat)
}
