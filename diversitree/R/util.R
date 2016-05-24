protect <- function(f, fail.value.default=NULL) {
  function(..., fail.value=fail.value.default, finite=TRUE) {
    if ( is.null(fail.value) )
      f(...)
    else {
      ret <- try(f(...), silent=TRUE)
      failed <- (inherits(ret, "try-error") ||
                 (finite && !is.finite(ret)))
      if ( failed )
        fail.value
      else
        ret
    }
  }
}

invert <- function(f) function(...) -f(...)

## Box constraints
boxconstrain <- function(f, lower, upper, fail.value=-Inf) {
  function(x, ...) {
    if ( any(x < lower | x > upper) )
      fail.value
    else
      f(x, ...)
  }
}

big.brother <- function(f, interval=1, file="") {
  f <- f # force argument to prevent recursion (pass by value)
  .x.eval <- list()
  .y.eval <- list()
  if ( interval < 0 )
    stop("Interval must be >= 0")
  ret <- function(x, ...) {
    i <- length(.x.eval) + 1
    if ( interval > 0 && i %% interval == 0 ) {
      cat(sprintf("[%s]", paste(formatC(x, 5), collapse=", ")),
          file=file)
      on.exit("\t -> [calculation failure]\n")
    } else if (interval > 0 ) {
      cat(".", file=file)
    }
    .x.eval[[i]] <<- x
    .y.eval[[i]] <<- ans <- f(x, ...)
    if ( interval > 0 && i %% interval == 0 ) {
      cat(sprintf("\t -> %2.5f\n", ans), file=file)
      on.exit()
    }
    ans
  }
  class(ret) <- c("big.brother", "function")
  attr(ret, "func") <- f
  ret
}

count.eval <- function(f) {
  n <- 0
  function(...) {
    n <<- n + 1
    f(...)
  }
}

set.defaults <- function(f, ..., defaults=NULL) {
  dots <- match.call(expand.dots=FALSE)[["..."]]
  if ( missing(defaults) )
    defaults <- dots
  else if ( is.list(defaults) )
    defaults <- c(dots, defaults)
  else
    stop("'defaults' must be a list")

  if ( is.null(defaults) )
    return(f)
  if ( !all(names(defaults) %in% names(formals(f))) )
    stop("Unknown defaults")
  att <- attributes(f)
  formals(f)[names(defaults)] <- defaults
  attributes(f) <- att[names(att) != "srcref"]
  f
}


quadratic.roots <- function(a, b, c)
  (-b + c(-1, 1) * sqrt(b*b - 4*a*c))/(2 * a)

discretize <- function(x, n, r=range(x)) {
  at <- seq(r[1], r[2], length.out=n+1)
  as.integer(cut(x, at, include.lowest=TRUE, labels=FALSE))
}

## Convert a matrix to a list by row.
matrix.to.list <- function(m) {
  n <- nrow(m)
  out <- vector("list", n)
  for ( i in seq_len(n) )
    out[[i]] <- m[i,]
  out
}

matrix.to.list <- function(m) {
  storage.mode(m) <- "double"
  .Call("matrix_to_list", m)
}

argnames.twopart <- function(base, n.level) {
  levels <- seq_len(n.level)
  paste(base, rep(levels, each=length(base)), sep=".")
}

## Construct a block diagonal matrix from a set of matrices in '...'.
## dimension names are assumed to be present or absent for all
## matrices without checking!
block.diagonal <- function(...) {
  matrices <- list(...)
  n <- sapply(matrices, dim)
  out <- matrix(0, sum(n[1,]), sum(n[2,]))
  rownames(out) <- unlist(lapply(matrices, rownames))
  colnames(out) <- unlist(lapply(matrices, colnames))
  offset <- c(0, 0)
  for ( i in seq_along(matrices) ) {
    out[seq(offset[1]+1, length.out=n[1,i]),
        seq(offset[2]+1, length.out=n[2,i])] <- matrices[[i]]
    offset <- offset + n[,i]
  }

  out
}

toC.int <- function(x) {
  x <- x - 1
  storage.mode(x) <- "integer"
  x
}

## Used by bm.vcv
dmvnorm2 <- function(x, mean, sigma, sigma.inv, log=FALSE) {
  distval <- mahalanobis(x, center=mean, cov=sigma.inv, TRUE)
  logdet <- as.numeric(determinant.matrix(sigma, TRUE)$modulus)
  logretval <- -(length(x) * log(2 * pi) + logdet + distval)/2
  if ( log )
    logretval
  else
    exp(logretval)
}

dmvnorm3 <- function(x, mean, sigma.inv, logdet, log=FALSE) {
  distval <- mahalanobis(x, center=mean, cov=sigma.inv, TRUE)
  logretval <- -(length(x) * log(2 * pi) + logdet + distval)/2
  if ( log )
    logretval
  else
    exp(logretval)
}

combine <- function(liks) {
  if ( !is.list(liks) )
    stop("liks must be a list")
  if ( !all(sapply(liks, inherits, "dtlik")) )
    stop("All elements of 'liks' must be diversitree likelihood functions")
  if ( length(unique(lapply(liks, class))) != 1 )
    stop("All functions must have the same class")
  if ( length(unique(lapply(liks, argnames))) != 1)
    stop("All functions must take the same argnames")
  if ( is.constrained(liks[[1]]) )
    stop("Cannot yet combine constrained functions -- constrain afterwards")
  ret <- function(pars, ...)
    sum(sapply(liks, function(f) f(pars, ...)))
  attributes(ret) <- attributes(liks[[1]])
  class(ret) <- c("combined", class(ret))
  info <- get.info(liks[[1]])
  info$phy <- NULL
  info$name.pretty <- sprintf("%s (combined [%d functions])",
                              info$name.pretty, length(liks))
  set.info(ret, info)
  ret
}

run.cached <- function(filename, expr, regenerate=FALSE) {
  if ( file.exists(filename) && !regenerate ) {
    readRDS(filename)
  } else {
    res <- eval.parent(substitute(expr))
    saveRDS(res, file=filename)
    res
  }
}
