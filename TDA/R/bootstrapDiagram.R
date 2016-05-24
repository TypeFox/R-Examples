bootstrapDiagram <- 
function(X, FUN, lim, by, maxdimension = length(lim) / 2 - 1,
         sublevel = TRUE, library = "Dionysus", B = 30, alpha = 0.05,
         distance = "bottleneck", dimension = min(1, maxdimension),
         p = 1, parallel = FALSE, printProgress = FALSE, weight = NULL, ...) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (class(FUN) != "function") {
    stop("FUN should be function")
  }
  if (!is.numeric(lim) || length(lim) %% 2 != 0) {
    stop("lim should be either a numeric matrix or a numeric vector of even elements")
  }
  if (!is.numeric(by) || any(by <= 0)) {
    stop("by should be positive")
  }
  if (2 * NCOL(X) != length(lim)) {
    stop("dimension of X does not match with lim")
  }
  if (length(by) != 1 && length(by) != NCOL(X)) {
    stop("by should be either a number or a vector of length equals dimension of grid")
  }
  if (!is.numeric(maxdimension) ||
      length(maxdimension) != 1 || maxdimension < 0) {
    stop("maxdimnsion should be a nonnegative integer")
  }
  if (!is.logical(sublevel)) {
    stop("sublevel should be logical")
  }
  if (library == "dionysus" || library == "DIONYSUS") {
    library = "Dionysus"
  }
  if (library == "phat" || library == "Phat") {
    library = "PHAT"
  }
  if (library != "Dionysus" && library != "PHAT") {
    stop("library should be a string: either 'Dionysus' or 'PHAT'")
  }
  if (!is.numeric(B) || length(B) != 1 || B < 1) {
    stop("B should be a positive integer")
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha > 1) {
    stop("alpha should be a number between 0 and 1")
  }
  if (distance != "wasserstein" && distance != "bottleneck") {
    stop("distance should be a string: either 'bottleneck' or 'wasserstein'")
  }
  if (!is.numeric(dimension) ||
      any(dimension < 0 || dimension > maxdimension)) {
    stop("dimension should be a integer or a vector of integer, with the range between 0 and maxdimension")
  }
  if (!is.numeric(p) || length(p) != 1 || p < 1) {
    stop("p should be a positive integer")
  }
  if (!is.logical(parallel)) {
    stop("parallel should be logical")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }
  if (!is.null(weight) && (!is.numeric(weight) ||
      (length(weight) != 1 && length(weight) != NROW(X)))) {
    stop("weight should be either NULL, a number, or a vector of length equals the number of sample")
  }

  X <- as.matrix(X)
  maxdimension <- min(maxdimension, NCOL(X) - 1)

  if (is.null(weight)) {
    Diag <- gridDiag(X = X, FUN = FUN, lim = lim, by = by, FUNvalues = NULL,
        maxdimension = maxdimension, sublevel = sublevel, library = library,
        location = FALSE, printProgress = FALSE, diagLimit = NULL, ...
      )[["diagram"]]
    if (distance == "wasserstein") {
      boostFUN <- function(i) {
        I <- sample(NROW(X), replace = TRUE, size = NROW(X))
        Diag1 <- gridDiag(X = X[I, , drop = FALSE], FUN = FUN, lim = lim,
            by = by, FUNvalues = NULL, maxdimension = maxdimension,
            sublevel = sublevel, library = library, location = FALSE,
            printProgress = FALSE, diagLimit = NULL, ...)[["diagram"]]
        width1 <- wasserstein(Diag, Diag1, p = p, dimension = dimension)
        if (printProgress) {
          cat(i, " ")
        }
        return(width1)
      }
    } else {
      boostFUN <- function(i) {
        I <- sample(NROW(X), replace = TRUE, size = NROW(X))
        Diag1 <- gridDiag(X = X[I, , drop = FALSE], FUN = FUN, lim = lim,
            by = by, FUNvalues = NULL, maxdimension = maxdimension,
            sublevel = sublevel, library = library, location = FALSE,
            printProgress = FALSE, diagLimit = NULL, ...)[["diagram"]]
        width1 <- bottleneck(Diag, Diag1, dimension = dimension)
        if (printProgress) {
          cat(i, " ")
        }
        return(width1)
      }
    }

  } else {
    Diag <- gridDiag(X = X, FUN = FUN, lim = lim, by = by, FUNvalues = NULL,
        maxdimension = maxdimension, sublevel = sublevel, library = library,
        location = FALSE, printProgress = FALSE, diagLimit = NULL,
        weight = weight, ...)[["diagram"]]
    if (distance == "wasserstein") {
      boostFUN <- function(i) {
        weightBoost <- rMultinom(size = sum(weight), prob = weight)
        Diag1 <- gridDiag(X = X, FUN = FUN, lim = lim, by = by,
            FUNvalues = NULL, maxdimension = maxdimension, sublevel = sublevel,
            library = library, location = FALSE, printProgress = FALSE,
            diagLimit = NULL, weight = weightBoost, ...)[["diagram"]]
        width1 <- wasserstein(Diag, Diag1, p = p, dimension = dimension)
        if (printProgress) {
          cat(i, " ")
        }
        return(width1)
      }
    } else {
      boostFUN <- function(i) {
        weightBoost <- rMultinom(size = sum(weight), prob = weight)
        Diag1 <- gridDiag(X = X, FUN = FUN, lim = lim, by = by, 
            FUNvalues = NULL, maxdimension = maxdimension, sublevel = sublevel,
            library = library, location = FALSE, printProgress = FALSE,
            diagLimit = NULL, weight = weightBoost, ...)[["diagram"]]
        width1 <- bottleneck(Diag, Diag1, dimension = dimension)
        if (printProgress) {
          cat(i, " ")
        }
        return(width1)
      }
    }
  }

  if (parallel) {
    boostLapply <- parallel::mclapply
  } else {
    boostLapply <- lapply
  }

  if (printProgress) {
    cat("Bootstrap: ")
  }
  width <- boostLapply(seq_len(B), FUN = boostFUN)
  if (printProgress) {
    cat("\n")
  }
  width <- stats::quantile(unlist(width), 1 - alpha)

  return (width)
}