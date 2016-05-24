# These functions are Copyright (C) 1998-2013 T. W. Yee   All rights reserved.

# Trying to get simulate() to work for some VGAM family functions.
# 20131228


# Last modified:
# 20131228: adapting simulate.vglm() from stats:::simulate.lm
# It comes from R 3.0.2.



simulate.vlm <- function (object, nsim = 1, seed = NULL, ...) {
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if (is.null(seed)) {
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  ftd <- fitted(object)
  nm <- names(ftd)
  n <- length(ftd)
  ntot <- n * nsim
  Fam <- if (inherits(object, "vlm")) {
#   object@family$family
    object@family
  } else {
#  "gaussian"
   stop("cannot get at the 'family' slot")
  }
#   if (!is.null(Fam@simslot)) {
#print("Hi1")
  val <-
    if (length(Fam@simslot) > 0) {
      Fam@simslot(object, nsim)
    } else {
      stop(gettextf("family '%s' not implemented", Fam), domain = NA)
    }
#print("val")
#print( val )
#stop("hello")
  if (!is.list(val)) {
    dim(val) <- c(n, nsim)
    val <- as.data.frame(val)
  } else {
    class(val) <- "data.frame"
  }
  names(val) <- paste("sim", seq_len(nsim), sep = "_")
  if (!is.null(nm)) 
    row.names(val) <- nm
  attr(val, "seed") <- RNGstate
  val
}





