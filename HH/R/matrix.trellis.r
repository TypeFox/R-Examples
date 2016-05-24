## turns trellis object x, with
##     length(dim(x))==1
## to a one-row trellis object xr with
##     length(dim(xr))==2 && dim(xr)[1] == dim(x)

"as.matrix.trellis" <- function(x, ..., row=FALSE, yname) {
  if (missing(yname))
    yname <- strsplit(deparse(x[[1]]), " ~ ")[[1]][1]
  xnames <- x$condlevels[[1]]
  nc <- length(xnames)

  x$index.cond <- list(1, 1:nc)
  x$condlevels <- list(yname, xnames)
  x$perm.cond  <- 1:2
  x$layout <- c(1, nc)

  if (row) {
  x$index.cond <- rev(x$index.cond)
  x$condlevels <- rev(x$condlevels)
  x$layout <- rev(x$layout)
  }
  x
}


## turns trellis object x, with
##     length(dim(x))==1
## to a one-row trellis object xr with
##     length(dim(xr))==2 && dim(xr)[1] == dim(x)
## matrix is not generic, hence this isn't really a method
## and hence I am using x, not data as the input object name



## turns trellis object x, with
##     length(dim(x))==1
## to a two-dimension trellis object xr

"matrix.trellis" <- function(x = NA, nrow = 1, ncol = 1,
                             byrow = FALSE, dimnames = NULL) {

  if (!byrow) {
    tmp <- ncol; ncol <- nrow; nrow <- tmp
  }
  if (byrow)
    dimnames <- rev(dimnames)

  x$index.cond <- list(1:ncol, 1:nrow)
  x$condlevels <-
      if (is.null(dimnames))
        x$index.cond
      else dimnames
  x$perm.cond  <- 1:2
  x$layout <- c(ncol, nrow)

  if (byrow) x else {
    warning("When byrow=FALSE,\n   this function produces a trellis object with the
right appearance,\n   but further processing with useOuterStrips and combineLimits will\n   a. ignore layout for non-square dim of trellis objects\n   b. give warning messages\n   c. display unreadable graphs.")
    ## recover()
    y <- x
    y$perm.cond <- rev(x$perm.cond)
    ## y$index.cond <- rev(x$index.cond)
    ## y$condlevels <- rev(x$condlevels)
    y$layout <- rev(x$layout)
    y
  }
}


## Take a two-dimensional trellis object and convert it to a
## one-dimensional trellis object
as.vector.trellis <- function(x, mode="any") {
  ## we are hijacking the variable 'mode' to use for the condlevels
  NN <- 1:prod(dim(x))
  x$index.cond <- list(NN)
  if (length(mode) != length(NN)) mode <- as.character(NN)
  x$condlevels <- list(mode)
  x$perm.cond  <- 1
  x
}


## Combine limits on a one-dimensional trellis object.  Do so by
## converting it two dimensions, combine the limits, and convert it
## back to one dimension,
combineLimits.trellisvector <-
  function(x, margin.x=1:2, margin.y=1:2,
           layout=x$layout,
           ncol=x$layout[1], nrow=x$layout[2],
           condlevels=x$condlevels[[1]],
           ...) {
    if (is.null(ncol)) ncol <- dim(x)
    if (is.null(nrow)) nrow <- 1
    byrow <- TRUE
    result <- as.vector.trellis(
      combineLimits(matrix.trellis(x, nrow=nrow, ncol=ncol, byrow=byrow),
                    margin.x= margin.x, margin.y= margin.y),
      condlevels)
    dimnames(result) <- dimnames(x)
    result
}



reorder.trellis <- function(x, X, ...) {
  y <- as.vector(x)
  y$panel.args <- y$panel.args[X]
  y
}
