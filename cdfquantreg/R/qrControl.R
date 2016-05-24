#' @title Control Optimization Parameters for CDF-Quantile Probability Distributions
#' @description Control Optimization Parameters for CDF-Quantile Probability Distributions.
#' @aliases cdfqr.control
#' @export
#' @param method Characters string specifying the method argument passed to \link[stats]{optim}.
#' @param maxit Integer specifying the maxit argument (maximal number of iterations) passed to \link[stats]{optim}.
#' @param hessian logical value, which indicates whether the numerical Hessian matrix from the optim output be used for estimation of the covariance matrix. Currently, only 'TRUE' is allowed, as only the analytical solution is employed.
#' @param trace Logical or integer controlling whether tracing information on the progress of the optimization should be produced
#' @return A list with the arguments specified.
#' @examples
#' 
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2', 't2', 
#' data = JurorData,control = cdfqr.control(trace = TRUE))


cdfqr.control<- function(method = "BFGS", maxit = 5000,
              hessian = TRUE, trace = FALSE){
  
  val <- list(method = method, maxit = maxit, hessian = hessian, 
               trace = trace)
  
  val
}

