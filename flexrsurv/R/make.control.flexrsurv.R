control.flexrsurv.additive <- function (epsilon = 1e-08, maxit = 25, trace = FALSE, ...){
  # controls for flexrsurv additive models are glm.controls
  return(glm.control(epsilon = epsilon, maxit = maxit, trace = trace))
  }



control.flexrsurv.multiplicative <- function (epsilon = 1e-8, maxit = 100, trace = FALSE, epsilon.glm = 1e-1, maxit.glm=30, ...){
  # controls for flexrsurv multiplicative models
  # epsilon : positive convergence tolerance of the iterative alternative algorithm;
  #           the iterations converge when |loglik - loglik_old}/(|loglik| + 0.1) < epsilon.
  # maxit  : integer giving the maximal number of iteration of the iterative alternative algorithm.
  # trace : logical indicating if output should be produced for each iteration.  
  # epsilon.glm,  maxit.glm  : contols passed to the inner glm call of the itertive alternative algorithm

    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("'epsilon' must be > 0", call.=FALSE)
    if (!is.numeric(epsilon.glm) || epsilon.glm <= 0) 
        stop("'epsilon.glm' must be > 0", call.=FALSE)
    if (!is.numeric(maxit.glm) || maxit.glm <= 0) 
        stop("maximum number of iterations 'maxit.glm' must be > 0", call.=FALSE)
    return(list(epsilon = epsilon, maxit = maxit, trace = trace, epsilon.glm = epsilon.glm, maxit.glm = maxit.glm))
  }



