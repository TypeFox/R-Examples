summary.bayescomm <-
function (object, chain, ...) {
  bsp <- paste("B", names(object$trace$B), sep = "$")
  if (chain == "R") {
    summary(mcmc(object$trace[[chain]], start = object$call$start,
                 thin = object$call$thin))
  } else if (chain %in% bsp) {
    summary(mcmc(object$trace$B[[substr(chain, 3, nchar(chain))]], start = object$call$start,
                 thin = object$call$thin))
  } else {
    stop("chain must be either 'R' or 'B$sp' for a named species sp")
  }
}