# ------------------------------------
# Author: Andreas Alfons
#         Erasmus University Rotterdam
# ------------------------------------

## add default column names to matrix
addColnames <- function(x) {
  # 'x' needs to be a matrix
  if(is.null(colnames(x))) colnames(x) <- paste("x", seq_len(ncol(x)), sep="")
  x
}

## add intercept column to design matrix
addIntercept <- function(x, check = FALSE) {
  if(!check || is.na(match("(Intercept)", colnames(x)))) {
    cbind("(Intercept)"=rep.int(1, nrow(x)), x)
  } else x
}

# ## backtransform regression coefficients to original scale (including intercept)
# backtransform <- function(beta, muY, sigmaY, mu, sigma) {
#   apply(beta, 2, 
#         function(b) {
#           b <- b * sigmaY / sigma
#           a <- muY - sum(b * mu)  # intercept
#           c("(Intercept)"=a, b)
#         })
# }

## call C++ back end
#' @useDynLib robustHD
callBackend <- function(..., PACKAGE) {
  # check the platfrom and if the RcppEigen back end is available
  # (RcppEigen back end does not work with 32-bit Windows)
  if(!isTRUE(.Platform$OS.type == "windows" && .Platform$r_arch == "i386") && 
       exists(".CallSparseLTSEigen")) {
    # RcppEigen back end from package sparseLTSEigen
    callFun <- get(".CallSparseLTSEigen")
    callFun(...)
  } else .Call(..., PACKAGE="robustHD")  # RcppArmadillo back end
}

## check the number of predictors to sequence for robust and groupwise LARS
## sequence predictors as long as there are twice as many observations
checkSMax <- function(sMax, n, p, groupwise = FALSE) {
  sMax <- rep(as.integer(sMax), length.out=1)
  if(groupwise) {
    m <- length(p)
    bound <- min(m, if(is.na(sMax)) floor(n/(2*mean(p))) else n-1)
  } else bound <- min(p, if(is.na(sMax)) floor(n/2) else n-1)
  if(!isTRUE(is.finite(sMax)) || sMax > bound) sMax <- bound
  sMax
}

## check the steps along the sequence for robust and groupwise LARS
checkSRange <- function(s, sMax = NA) {
  s <- as.integer(s)
  if(length(s) == 0) s <- c(0, sMax)
  else if(length(s) == 1) s <- rep.int(s, 2)
  else s <- s[1:2]
  if(!isTRUE(is.finite(s[1]))) s[1] <- 0
  if(isTRUE(is.finite(s[2]))) {
    if(s[1] > s[2]) s[1] <- s[2]
  } else s[2] <- NA
  s
}

## check steps for coef(), fitted(), residuals(), predict(), ... methods
checkSteps <- function(s, sMin, sMax, recycle = FALSE, ...) {
  ok <- is.numeric(s) && length(s) > 0 && all(is.finite(s))
  if(ok) {
    if(isTRUE(recycle)) {
      s[s < sMin] <- sMin
      s[s > sMax] <- sMax
    } else ok <- all(s >= sMin) && all(s <= sMax)
  }
  if(!ok) stop(sprintf("invalid step, must be between %d and %d", sMin, sMax))
  s
}

## copy names from a vector or matrix to another vector or matrix
copyNames <- function(from, to, which = "col", target = "row") {
  # read names from source
  if(is.null(dim(from))) nam <- names(from) 
  else if(which == "row") nam <- rownames(from)
  else if(which == "col") nam <- colnames(from)
  # write names to target
  if(is.null(dim(to))) names(to) <- nam
  else if(target == "row") rownames(to) <- nam
  else if(target == "col") colnames(to) <- nam
  # return object
  to
}

## utility function to get default labels for plot
defaultLabels <- function(x) UseMethod("defaultLabels")

defaultLabels.seqModel <- defaultLabels.sparseLTS <- function(x) {
  as.character(seq_along(removeIntercept(coef(x))))
}

defaultLabels.grplars <- function(x) {
  assign <- x$assign
  labels <- split(as.character(assign), assign)
  p <- sapply(labels, length)  # number of variables per group
  append <- which(p > 1)
  labels[append] <- mapply(function(l, p) paste(l, seq_len(p), sep="."), 
                           labels[append], p[append], SIMPLIFY=FALSE)
  unsplit(labels, assign)
}

## utility function to get default main plot title
defaultMain <- function() "Coefficient path" 

## drop dimension in case of matrix with one column
dropCol <- function(x) {
  d <- dim(x)
  if(is.null(d[2]) || d[2] != 1) x
  else if(d[1] == 1) {
    # drop() drops all names for a 1x1 matrix
    names <- rownames(x)
    x <- drop(x)
    names(x) <- names
    x
  } else drop(x)
}

## construct blocks of original and lagged values for time series models
fitBlocks <- function(x, y, h = 1, p = 2, intercept = FALSE) {
  n <- length(y)
  tsBlocks(x, y, p=p, subset=-((n-h+1):n), intercept=intercept)
}

# ## find argument names of functions
# findArgNames <- function(..., removeDots = TRUE) {
#   argNames <- unique(unlist(lapply(list(...), function(f) names(formals(f)))))
#   if(removeDots) {
#     argNames <- setdiff(argNames, "...")
#   }
#   argNames
# }

## find indices of h smallest observations
findSmallest <- function(x, h) {
  # call C++ function
  callBackend("R_findSmallest", R_x=as.numeric(x), R_h=as.integer(h))
}

## compute coefficients of hyperplane through data points
hyperplane <- function(x) {
  p <- ncol(x)
  y <- -x[, p]  # right hand side
  x <- x[, -p, drop=FALSE]
  tx <- t(x)
  theta <- solve(tx %*% x) %*% tx %*% y
  c(theta, 1)
}

## obtain the degrees of freedom of a model (number of nonzero parameters)
modelDf <- function(beta, tol = .Machine$double.eps^0.5) {
  length(which(abs(beta) > tol))
}

## construct blocks of original and lagged values for prediction from time 
## series models
newdataBlocks <- function(x, y, h = 1, p = 2, intercept = TRUE) {
  n <- length(y)
  tsBlocks(x, y, p=p, subset=(n-h-p+2):n, intercept=intercept)
}

## find indices of h smallest observations
partialOrder <- function(x, h) {
  # call C++ function
  callBackend("R_partialOrder", R_x=as.numeric(x), R_h=as.integer(h))
}

# ## find indices of h smallest observations
# partialSort <- function(x, h) {
#   # call C++ function
#   callBackend("R_partialSort", R_x=as.numeric(x), R_h=as.integer(h))
# }

## remove intercept column from design matrix
removeIntercept <- function(x, pos) {
  haveVector <- is.null(dim(x))
  if(missing(pos)) {
    names <- if(haveVector) names(x) else colnames(x)
    pos <- match("(Intercept)", names, nomatch=0)
  }
  if(pos > 0) {
    if(haveVector) x[-pos] else x[, -pos, drop=FALSE]
  } else x
}

# ## prepend something to column names of a matrix
# prependColnames <- function(x, prefix) {
#   colnames(x) <- paste(prefix, colnames(x), sep=".")
#   x
# }
