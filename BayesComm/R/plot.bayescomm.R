plot.bayescomm <-
function (x, chain, ...) {
  bsp <- paste("B", names(x$trace$B), sep = "$")
  if (chain =="R") {
    plot(mcmc(x$trace[[chain]], start = x$call$start,
              thin = x$call$thin))
  } else if (chain %in% bsp) {
    plot(mcmc(x$trace$B[[substr(chain, 3, nchar(chain))]],
              start = x$call$start, thin = x$call$thin))
  } else {
    stop("chain must be either 'R' or 'B$sp' for a named species sp")
  }
}