Mort2Dsmooth_optimize <-
function(x, y, Z, offset,
                                  wei, psi2,
                                  Bx, By,
                                  nbx, nby,
                                  RTBx, RTBy,
                                  Px, Py,
                                  a.init,
                                  MON, TOL1, TOL2,
                                  RANGEx, RANGEy,
                                  MAX.IT, MET){
  ## the function first applies a rough grid search
  ## then around the optimal values,
  ## it search with cleversearch
  ## a more precise value for the smoothing parameters

  ## Input:
  ## x: abcissa of data
  ## y: ordinate of data
  ## Z: matrix of count response
  ## offset: an a priori known component (optional)
  ## wei: a matrix of weights
  ##      to be used in the fitting process
  ## psi2: overdispersion parameter
  ## Bx: B-splines basis over x
  ## By: B-splines basis over y
  
  ## DtD: inner product of the matrix of differences
  ## a.init: starting coefficients
  ## MON: logical on monitoring
  ## TOL1: relative convergence tolerance
  ## TOL2: difference between two adjacent
  ##       smoothing parameters in the second grid search,
  ##       log-scale.
  ## RANGEx: range in between grid-search for x
  ## RANGEy: range in between grid-search for y
  ## MAX.IT: the maximum number of iterations
  ## MET: criterion for optimizing lambda (bic or aic)
  
  ## Output: a list containing
  ## lambdas: the optimal smoothing parameters

  ## object function 
  Mort2Dsmooth_opt_ic <- function(X){
    FIT <- Mort2Dsmooth_estimate(x=x, y=y, Z=Z,
                                 offset=offset,
                                 psi2=psi2, wei=wei,
                                 Bx=Bx, By=By,
                                 nbx=nbx, nby=nby,
                                 RTBx=RTBx, RTBy=RTBy,
                                 lambdas=c(X[1], X[2]),
                                 Px=Px, Py=Py,
                                 a.init=a.init,
                                 MON=MON,
                                 TOL1=TOL1, 
                                 MAX.IT=MAX.IT)
    return(ifelse(MET==2, FIT$aic, FIT$bic))
  }

  ## ROUGH SEARCHING
  ## possible lambdas
  lambdas.x0 <- seq(log10(RANGEx[1]),log10(RANGEx[2]),
                    TOL2*4)
  lambdas.y0 <-  seq(log10(RANGEx[1]),log10(RANGEx[2]),
                     TOL2*4)
  ## length of possible lambdas
  by.lambdas.x0 <- length(lambdas.x0)
  by.lambdas.y0 <- length(lambdas.y0)
  by.lambda.0 <- max(by.lambdas.x0, by.lambdas.y0)
  ## starting values for lambdas in cleversearch
  l.st.x0 <- median(lambdas.x0)
  l.st.y0 <- median(lambdas.y0)
  ## optimizing lambda
  lambdas.hat0 <- cleversearch(fn=Mort2Dsmooth_opt_ic,
                               lower=c(log10(RANGEx[1]),
                                       log10(RANGEy[1])),
                               upper=c(log10(RANGEx[2]),
                                       log10(RANGEy[2])),
                              startvalue=c(l.st.x0,l.st.y0),
                               ngrid=by.lambda.0,
                               logscale=TRUE,
                               verbose=FALSE)[[1]]
  
  ## FINE SEARCHING
  ## sub-optimal lambdas=starting values in cleversearch
  l.st.x <- log10(lambdas.hat0[1])
  l.st.y <- log10(lambdas.hat0[2])
  ## range in which search
  min.l.x <- max(l.st.x-TOL2*4, log10(RANGEx[1]))
  max.l.x <- min(l.st.x+TOL2*4, log10(RANGEx[2]))
  min.l.y <- max(l.st.y-TOL2*4, log10(RANGEy[1]))
  max.l.y <- min(l.st.y+TOL2*4, log10(RANGEy[2]))
  ## possible lambdas
  lambdas.x <- seq(min.l.x, max.l.x, TOL2)
  lambdas.y <- seq(min.l.y, max.l.y, TOL2)
  ## length of possible lambdas
  by.lambdas.x <- length(lambdas.x)
  by.lambdas.y <- length(lambdas.y)
  by.lambda <- max(by.lambdas.x, by.lambdas.y)
  ## optimizing lambda
  lambdas.hat <- cleversearch(fn=Mort2Dsmooth_opt_ic,
                              lower=c(min.l.x, min.l.y),
                              upper=c(max.l.x, max.l.y),
                              startvalue=c(l.st.x, l.st.y),
                              ngrid=by.lambda,
                              logscale=TRUE,
                              verbose=FALSE)[[1]]
  
  return(lambdas.hat)
}
