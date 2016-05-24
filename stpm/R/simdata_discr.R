#' Multi-dimension simulation function
#'@references Akushevich I., Kulminski A. and Manton K. (2005), Life tables with covariates: Dynamic model 
#'for Nonlinear Analysis of Longitudinal Data. Mathematical Population Studies, 12(2), pp.: 51-80.
#'<DOI:10.1080/08898480590932296>.
#' @param N Number of individuals
#' @param a A k by k matrix, which characterize the rate of the adaptive response.
#' @param f1 A particular state, which is a deviation from the normal (or optimal). This is a vector with length of k.
#' @param Q A matrix k by k, which is a non-negative-definite symmetric matrix.
#' @param f A vector-function (with length k) of the normal (or optimal) state.
#' @param b A diffusion coefficient, k by k matrix.
#' @param mu0 mortality at start period of time.
#' @param theta A displacement coefficient of the Gompertz function.
#' @param ystart A vector with length equal to number of dimensions used, defines starting values of covariates.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @param dt A time step (1 by default).
#' @return A table with simulated data.
#' @examples
#' library(stpm)
#' data <- simdata_discr(N=100, ystart=80)
#' head(data)
#'
simdata_discr <- function(N=100, a=-0.05, f1=80, Q=2e-8, f=80, b=5, mu0=1e-5, theta=0.08, ystart=80, tstart=30, tend=105, dt=1) {
  
  k <- length(ystart)
  
  if ( (dim(as.data.frame(a))[1] != k) & (dim(as.data.frame(a))[2] != k) ) {
    stop("Dimenstions if \'a\' are not equal to k. It must be a k x k matrix.")
  } else if( (dim(as.data.frame(Q))[1] != k) & (dim(as.data.frame(Q))[2] != k) ) {
    stop("Dimenstions of \'Q\' are not equal to k. It must be a k x k matrix.")
  } else if((dim(as.data.frame(f1))[2] != k)) { 
    stop("\'f1'\ must be a 1 x k matrix")
  } else if( (dim(as.data.frame(f1))[2] != k)){ 
    stop("\'f'\ must be a 1 x k matrix")
  } else if(dim(as.data.frame(b))[1] != k) {
    stop("\'b'\ must be a 1 x k matrix")
  } else if((dim(as.data.frame(ystart))[1] != k) ) {
    stop("\'y'\ must be a 1 x k vector")
  }
  
  # Re-calculating parameters:
  u_ <- matrix((f1 %*% (-1*a)), nrow=k, ncol=1)
  R_ <- matrix((diag(k) + a), nrow=k, ncol=k)
  Sigma_ <- matrix(b, nrow=k, ncol=1)
  mu0_ <- mu0 + f %*% Q %*% t(f)
  b_ <- matrix((-2*f %*% Q), nrow=k, ncol=1)
  Q_ <- matrix(Q, nrow=k, ncol=k)
  theta_ <- theta
  ystart = matrix(ystart, nrow=k, ncol=1)
  simulated = .Call("simdata_ND", N, u_, R_, Sigma_, mu0_, b_, Q_, theta_, tstart, ystart, tend, k, dt);
  
  data_names <- c()
  for(n in 1:k) {
    data_names <- c(data_names, paste("y",n, sep=''), paste("y",n, ".next", sep=''))
  }
  colnames(simulated) <- c("id", "xi", "t1", "t2", data_names)
  
  invisible(simulated)
}

