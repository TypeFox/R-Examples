
rvnsims <- function (x) {
  UseMethod("rvnsims")
}

rvnsims.rv <- function (x) {
  sapply(unclass(x), length)
}

rvnsims.rvsummary <- function (x) {
  unlist(rvattr(x, "n.sims"), use.names=TRUE)
}

rvnsims.default <- function (x) {
  if (!(is.atomic(x) || is.recursive(x))) {
    stop("rvnsims: no applicable method for class '", class(x), "'")
  }
  rep.int(1, length(x))
}

setnsims <- function (n.sims) {
  ## setnsims - get or set the default number of simulations (a global variable)
  if (length(n.sims)>0 && is.numeric(n.sims) && (!is.na(n.sims[1])) && n.sims[1]>=2) {
    n.sims <- as.integer(ceiling(n.sims[1]))
    oldn.sims <- rvpar("n.sims")
    rvpar(n.sims=n.sims)
  } else {
    stop('Invalid number of simulations (must be at least 2)', n.sims[1])
  }
  return(oldn.sims)
}

getnsims <- function () {
  n.sims <- rvpar("n.sims")
  if (!is.integer(n.sims) || n.sims<2) {
    stop("Invalid number of simulations - please set with setnsims(...)")
  }
  return(n.sims)
}

