.margins = function(x, margins){
  if (length(margins) == 1) {
    x = apply(x, 2, .one.mar, margins = margins)
  }
  else {
    for (i in 1:NCOL(x)) {
      x[, i] = .one.mar(x[, i], margins[i])
    }
  }
  x
}

.one.mar = function(x, margins){
  n = NROW(x)
  if (margins == "ranks") {
    ecdf(x)(x)*n/(n+1)
  }
  else {
    boundary = 10000
    if ((margins == "beta") | (margins == "cauchy") | (margins == "chisq") | 
          (margins == "f") | (margins == "gamma") | (margins == "lnorm") | 
          (margins == "norm") | (margins == "t") | (margins == "weibull")) {
      loglik = function(par, x) {
        sum(log(eval(do.call(paste("d", margins, sep = ""), 
                             args = list(x = x, par[1], par[2])))))
      }
      op = constrOptim(theta = c(1, 1), f = loglik, grad = NULL, 
                       ui = matrix(c(1, 0, -1, 0, 0, 1), nrow = 3, byrow = TRUE), 
                       ci = c(-rep(boundary, 2), 0), x = x, control = list(fnscale = -1), 
                       hessian = FALSE)
      eval(do.call(paste("p", margins, sep = ""), args = list(q = x, 
                                                           op$par[1], op$par[2])))
    }
    else {
      if ((margins == "exp")) {
        op = optimise(f = function(par, x) {
          sum(log(eval(do.call(paste("d", margins, sep = ""), 
                               args = list(x = x, par)))))
        }, x = x, lower = 1e-04, upper = 100, maximum = TRUE)$maximum
        eval(do.call(paste("p", margins, sep = ""), args = list(q = x, 
                                                             op)))
      }
    }
  }
}