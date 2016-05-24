powell.control <- function(trace = 0, rhobeg = log(100)/10, rhoend = 1e-4,
                           maxit = 1e5, fnscale = 1., parscale = 1., ...) {
  list(trace = as.integer(trace), rhobeg = rhobeg, rhoend = rhoend,
       maxit = maxit, fnscale = fnscale, parscale = parscale)
}
