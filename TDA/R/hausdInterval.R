hausdInterval <- 
function(X, m, B = 30, alpha = 0.05, parallel = FALSE, printProgress = FALSE) {
     
  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(m) || length(m) != 1 || m < 0) {
    stop("m should be a nonnegative integer")
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    stop("B should be a positive integer")
  }
  if (!is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha should be a number between 0 and 1")
  }
  if (!is.logical(parallel)) {
    stop("parallel should be logical")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }

  X <- as.matrix(X)
  n <- nrow(X)
  width <- rep(0, B)
  if (parallel) {
    boostLapply <- parallel::mclapply
  } else {
    boostLapply <- lapply
  }

  if (printProgress) {
    cat("Subsampling: ")
  }
  width <- boostLapply(1:B, FUN = function(i) {
      I <- sample(1:n, replace = FALSE, size = m)
      Y <- as.matrix(X[I,])
      LL <- max(FNN::knnx.dist(Y, X, k = 1, algorithm = "kd_tree"))
      if (printProgress) {
        cat(i," ")
      }
      return(LL)
    })
  if (printProgress) {
    cat("\n")
  }
  width <- unlist(width)
  width <- 2 * stats::quantile(width, 1 - alpha)

  out <- width

  return(out)
}