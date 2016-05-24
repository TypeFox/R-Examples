#' @importFrom stats deviance df.residual fitted rnorm runif
#' @keywords internal
simulate.nls <- function (object, nsim = 1, seed = NULL, ...) {
  
  # Initialize random number generator
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  if (is.null(seed)) 
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  
  # Simulate new response values
  ftd <- fitted(object)  # fitted values
  nm <- names(ftd)  # row names
  n <- length(ftd)  # number of observations
  ntot <- n * nsim  # number of observations to simulate
  vars <- deviance(object) / df.residual(object)  # residual variance
  if (!is.null(object$weights)) {
    vars <- vars/object$weights
  }
  val <- ftd + rnorm(ntot, sd = sqrt(vars))  # simulated response values
  
  # Return simulated response values with appropriate dimension, class, etc.
  if (!is.list(val)) {
    dim(val) <- c(n, nsim)  # assign dimensions (should b n-by-nsim)
    val <- as.data.frame(val)  # convert to data frame
  } else {
    class(val) <- "data.frame"
  }
  names(val) <- paste("sim", seq_len(nsim), sep = "_")
  if (!is.null(nm)) {
    row.names(val) <- nm
  }
  attr(val, "seed") <- RNGstate
  val
  
}
