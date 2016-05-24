#General utilities - not exported

#' @import stats
#' @import utils

#Strange errors were reported on solaris-sparc, this attempts to avoid them
safe.solve <- function(a, b) {
  if (missing(b)) {
    try.result <- try(x <- solve(a))
  } else {
    try.result <- try(x <- solve(a, b))
  }
  if (inherits(try.result, "try-error")) {
    if (missing(b)) {
      x <- matrix(nrow = nrow(a), ncol = ncol(a))
    } else {
      x <- matrix(nrow = ncol(a), ncol = ncol(b))
    }
    warning("Error in solve(), standard errors not available")
  }
  return(x)
}

#like seq, but returns integer(0) if from > to   (always increments by 1)
sseq <- function(from, to) {
  if (from > to) return(integer(0))
  seq(from, to)
}

#source: http://stackoverflow.com/questions/23274170/how-to-efficiently-check-if-a-matrix-is-in-binary-form-e-g-all-1s-or-0s
IsBinary <- function(mat) {
  identical(mat, as.numeric(as.logical(mat)))
}

# scale to 0.01, 0.99 and take logit
LogitScale <- function(x) {
  qlogis(Scale(x, 0.01, 0.99))
}

Scale <- function(x, min.y, max.y) {
  if (all(is.na(x))) stop("all NA in Scale")
  r <- range(x, na.rm = TRUE)
  if (diff(r) > 0) {
    return((x - r[1])/diff(r) * (max.y - min.y) + min.y)
  } else {
    #only one value of x
    if (r[1] >= min.y && r[1] <= max.y) {
      #if the one value is in [min.y, max.y], return it
      return(rep(r[1], length(x)))
    } else {
      #otherwise return mean(min.y, max.y)
      return(rep(mean(c(min.y, max.y)), length(x)))
    }
  }
}

# If x is a matrix, keep it; if x is a vector, make it a 1 column matrix
AsMatrix <- function(x) {
  if (is.matrix(x)) {
    return(x)
  } else if (is.vector(x)) {
    dim(x) <- c(length(x), 1)
    return(x)
  } else {
    stop("AsMatrix input should be a matrix or vector")
  }
}

is.equal <- function(...) {
  isTRUE(all.equal(...))
}

#same as model.matrix, but leave NA rows
model.matrix.NA <- function(object, data) {
  MM <- model.matrix(object, data)
  MM <- MM[match(rownames(data), rownames(MM)), , drop=FALSE]
  rownames(MM) <- rownames(data)
  return(MM)
}

#if warning is in ignoreWarningList, ignore it; otherwise post it as usual
SuppressGivenWarnings <- function(expr, warningsToIgnore) {
  h <- function (w) {
    if (w$message %in% warningsToIgnore) invokeRestart( "muffleWarning" )
  }
  withCallingHandlers(expr, warning = h )
}

get.stack <- function(x, ifnotfound=NULL, mode="any") {
  #Look for object x (string) searching through the calling stack, starting at the current parent frame and going back through the parents
  if (!is.character(x)) stop("x must be a character string")
  for (f in (sys.parent(1)):0) {
    if (exists(x, envir=sys.frame(f), mode=mode, inherits=FALSE)) {
      return(get(x, envir=sys.frame(f), mode=mode, inherits=FALSE))
    }
  }
  return(ifnotfound)
}


rexpit <- function(x) rbinom(n=length(x), size=1, prob=plogis(x))

# Given row and column numbers of a matrix with num.rows rows, compute the single index
sub2ind <- function(row, col, num.rows) {
  return((col - 1) * num.rows + row)
}

repmat <- function(X,m,n){
  #R equivalent of repmat (matlab)
  mx <- dim(X)[1]
  nx <- dim(X)[2]
  if ((m == 0) || (n == 0)) return(matrix(numeric(0), nrow=mx*m, ncol=nx*n)) #avoids warnings when m or n is 0
  return(matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T))
}

# from Ken Williams on StackOverflow
Mode <- function(x, na.rm=FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

drop3 <- function(x) {
  #if x is an array with 3 dimensions and third dimension has one level, return a matrix with it dropped; otherwise error
  return(dropn(x, 3))
}

dropn <- function(x, n) {
  #if x is an array with n dimensions and nth dimension has one level, return a matrix with it dropped; otherwise error
  stopifnot(length(dim(x))==n)
  stopifnot(dim(x)[n]==1)
  dn <- dimnames(x)
  dim(x) <- dim(x)[1:(n-1)]
  dimnames(x) <- dn[1:(n-1)]
  return(x)
}
