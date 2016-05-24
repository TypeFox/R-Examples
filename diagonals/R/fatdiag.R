#' Fat Matrix Diagonals
#'
#' @aliases fatdiag<-
#' @param x a matrix where the dimensions are integer multiples of size or integer dividors of steps
#' @param steps the required number of steps (block matrices) across the diagonal
#' @param size the width or height of the matrix being dropped over the diagonal of matrix x
#' @param nrow the number of rows
#' @param ncol the number of columns
#' @details Either steps or size is expected to be provided.
#' @export
#' @examples
#' fatdiag(12, steps=3)
#'
#' ( m <- matrix(111, nrow=6, ncol=9) )
#' fatdiag(m, steps=3) <- 5
#'
#' fatdiag(m, steps=3)
#'
#' fatdiag(12, size=4)
#'
#' fatdiag(12, size=c(3,4) )

fatdiag <- function( x = 1, steps=NULL, size=NULL, nrow=NULL, ncol=NULL) {

  if (length(size) == 1)
    size <- c(size, size)

  if (length(x) == 1) {

    if ( !is.null(nrow) && !is.null(ncol) ) {
      dx <- as.vector( c(nrow,ncol) )
    } else if ( !is.null(nrow) && is.null(ncol)) {
      if ("common denominator x and nrow")
        stop("nrow and x do not have a common denominator")
      dx <- as.vector(c(nrow, x))
    } else if ( is.null(nrow) && !is.null(ncol)) {
      if ("common denominator x and ncol")
        stop("ncol and x do not have a common denominator")
      dx <- c(x, ncol)
    } else if ( is.null(nrow) && is.null(ncol) && is.null(steps) ) {
      steps <- x %/% max(size)
      dx    <- size * steps
    } else {
      dx <- c(x, x)
    }

    if ( !all(dx %% steps == 0) )
      stop("steps is not an integer divisor of x on all dimensions")

    # create a fat diagonal matrix
    m <- matrix(0, nrow=dx[1], ncol = dx[2])
    fatdiag(m, steps = steps, size = size) <- 1
    return(m)

  } else if ( is.null( dim(x) ) && length(x > 1) ) {

    if ( !is.null(nrow) && !is.null(ncol) ) {
      dx <- as.vector( c(nrow,ncol) )
    } else if ( !is.null(nrow) && is.null(ncol)) {
      if ("common denominator x and nrow")
        stop("nrow and x do not havea common denominator")
      dx <- as.vector(c(nrow, x))
    } else if ( is.null(nrow) && !is.null(ncol)) {
      if ("common denominator x and ncol")
        stop("ncol and x do not have a common denominator")
      dx <- c(x, ncol)
    } else if ( is.null(nrow) && is.null(ncol) && is.null(steps) ) {
      steps <- length(x) %/% ( max(size) * min(size) )
      dx    <- size * steps
    } else {
      size <- c( sqrt(length(x)/steps), sqrt(length(x)/steps) )
      dx   <- c( (length(x) / size),  (length(x) / size)  )
    }

    if ( !all(dx %% steps == 0) )
      stop("steps is not an integer divisor of x on all dimensions")


    # create a fat diagonal matrix
    m <- matrix(0, nrow=dx[1], ncol = dx[2])
    fatdiag(m, steps = steps, size = size) <- x
    return(m)

  } else if ( length(dim(x)) == 2) {

    # extract the fat diagonal

    dx <- dim(x)
    size <- dx %/% steps
    # split dimension according to steps
    spl1 <- split_vector(1:dx[1], steps = steps)
    spl2 <- split_vector(1:dx[2], steps = steps)
    # create vectors
    a <- vector()
    b <- vector()
    for (i in 1:steps) {
      a <- c(a, rep(spl1[[i]], times = size[2]) )
      b <- c(b, rep(spl2[[i]], each  = size[1]) )
    }
    return( x[cbind(a,b)] )
  } else {
    stop("x is not a valid vector or matrix")
  }

}

#' @describeIn fatdiag the set version of fatdiag
#' @title fatdiag set
#' @aliases fatdiag
#' @param on_diagonal should the operation be apply to the elements on the fat diagonal.
#' @param value replacement value
#' @export
`fatdiag<-` <- function( x, steps = NULL, size = NULL, on_diagonal=TRUE, value ) {


  if (length(size) == 1)
    size <- c(size, size)

  # save dimensions
  dx <- dim(x)
  # save value length
  lv <- length(value)

  # square if dimensions are right
  if (length(dx) != 2L)
    stop("not a matrix")

  # determine the size of the step
  if ( is.null(steps) && !is.null(size) ) {

    # coerce to integer
    # size <- floor(size)

    # create steps
    steps <- max(dx) %/% max(size)

  } else if ( !is.null(steps) & is.null(size) ) {

    # calculate size
    size <- dx %/% steps

  } else if (is.null(steps) & is.null(size) ) {

    # issue warning
    warning("Both steps and size parameters are NULL, trying to guess size")

    if ( all(sqrt(dx) %% 1 == 0) ) {
      size  <- sqrt(dx)
      steps <- sqrt(dx[1])
      warning("using the square root as steps and size")
    } else {
      # set to unit
      size <- 1L
      # create steps
      steps <- dx[1] %/% size
    }

  }

  # check that dimensions of this matrix are a multiple of size
  if( !all(dx %% size == 0) )
    stop("Matrix dimensions are not a multiple of size")

  if (as.integer(size[1]*size[2]*steps) %% lv != 0 && lv != 1L)
    stop("value fat diagonal has wrong length")


  # split dimension according to steps
  spl1 <- split_vector(1:dx[1], steps = steps)
  spl2 <- split_vector(1:dx[2], steps = steps)

  # create vectors
  a <- vector()
  b <- vector()

  if (!on_diagonal) {
    for (i in 1:steps) {
      a <- c(a, rep(spl1[[i]], times=dx[1]-size) )
      b <- c(b, rep(setdiff(1:dx[1], spl2[[i]]), each=size) )
      }
    } else {
    # insert combinations
    for (i in 1:steps) {
      a <- c(a, rep(spl1[[i]], times = size[2]) )
      b <- c(b, rep(spl2[[i]], each  = size[1]) )
      }
    }

  # replace
  x[cbind(a,b)] <- value

  # return output
  return(x)

}
