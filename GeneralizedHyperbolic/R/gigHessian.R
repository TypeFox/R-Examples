### Calculate Hessian given fitted parameters
gigHessian <- function(x, param, hessianMethod = c("tsHessian", "exact"),
                       whichParam = 1) {
  if (hessianMethod == "exact") {
    stop("Exact hessian not available yet. Use method tsHessian instead.")
  }
  if (whichParam != 1){
    stop("Only parameterization 1 has been implemented at present")
  }
  if (hessianMethod == "tsHessian") {
    llfuncH <- function(param) {
      loggigDens <- param[3]/2*log(exp(param[2])/exp(param[1])) -
                    log(2*besselK(sqrt(exp(param[1])*exp(param[2])),
                                  nu = param[3])) + (param[3] - 1)*log(x) -
                    1/2*(exp(param[1])*x^-1 + exp(param[2])*x)
      as.numeric(loggigDens)
      return(sum(loggigDens))
    }
  }

  hessian <- tsHessian(param = param, fun = llfuncH)
  return(hessian)
}
