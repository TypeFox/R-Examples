.missingFormals <- function(f, g, offset=0) {
  fArgs <- formals(f)
  hasCatchAll <- !is.null(fArgs[["..."]])
  fArgs["..."] <- NULL
  skipIndices <- rep(FALSE, length(fArgs))
  if (offset > 0) {
    skipIndices[1:offset] <- TRUE
  }
  fArgs[skipIndices] <- NULL
  if (length(fArgs) == 0) {
    if (hasCatchAll) {
      return(NULL)
    } else {
      return(names(g))
    }
  } else {
    hasMissingFormals <- !(names(fArgs) %in% names(g))
    if (any(hasMissingFormals)) {
      return(names(fArgs[hasMissingFormals]))
    } else {
      return(NULL)
    }
  }
}

hermite <- function(x, n, prob=TRUE) {
  isBadDegree <- any(n < 0 || !isInteger(n))
  if(isBadDegree) {
    stop("Argument 'n' must be a vector of non-negative integers!")
  }
  isBadLength <- (length(n) != 1) && (length(x) != length(n)) && (length(x) != 1)
  if(isBadLength) {
    stop(paste("Argument 'n' must be either a vector of same length",
               "as argument 'x',\n  a single integer or 'x' must be a ",
               "single value!", sep=""))
  }
  H <- function(x, n) {
    if (n <= 1) {
      return(switch(n + 1, 1, x))
    } else {
      return(x * Recall(x, n - 1) - (n - 1) * Recall(x, n - 2))
    }
  }
  scale <- 1
  if (!prob) {
    x <- sqrt(2) * x
    scale <- 2 ^ (n / 2)
  }
  scale * mapply(H, x, n)
}
