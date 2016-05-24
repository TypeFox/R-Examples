
##*****************************************************************************
## Grid linear interpolation in arbitrary dimension through Cardinal
## Basis.
##
## This is a grid interpolation as in \code{\link{gridInt}} but it is
## required here that the one-dimensional interpolation method is
## \emph{linear} w.r.t. the vector of interpolated values.  For each
## dimension, a one-dimensional interpolation is carried out,
## leading to a collection of interpolation problems each with a
## dimension reduced by one. The \emph{same cardinal basis} can be
## used for all interpolations w.r.t. the same variable, so the
## number of call to the \code{interpCB} function is equal to the
## interpolation dimension, making this method \emph{very fast} compared to
## the general grid interpolation as implemented in
## \code{\link{gridInt}}, see the \bold{Examples} section.
## 
## @title Grid interpolation in arbitrary dimension through Cardinal
## Basis
##
## @param X An object that can be coerced into \code{Grid}. This
## can be a data.frame or a matrix in \emph{Scattered Data style} in
## which case the column number is equal to the spatial dimension
## \eqn{d}, and the row number is then equal to the number of nodes
## \eqn{n}. But it can also be a \code{Grid} object previously
## created.  A data frame or matrix \code{X} will first be coerced
## into \code{Grid} by using the the S3 method
## \code{\link{as.Grid}}.
## 
## @param Y Response to be interpolated. It must be a vector of
## length \eqn{n} equal to the number of nodes. When \code{X} has
## class \code{Grid}, the order of the elements in \code{Y} must
## conform to the order of the nodes as given in \code{X}, see the
## help for \code{\link{gridInt}}.
##
## @param Xout Interpolation locations. Can be a vector or a
## matrix. In the first case, the length of the vector must be equal
## to the spatial dimension \eqn{d} as given by \code{xLev} or
## \code{X}.  In the second case, each row will be considered as a
## response to be interpolated, and the number of columns of
## \code{Xout} must be equal to the spatial dimension.
##
## @param interpCB Function evaluating the interpolation Cardinal
## Basis. This function must have as its first 2 formals 'x', and
## 'xout'.  It must return a matrix with \code{length(x)} columns
## and \code{length(xout)} rows. The \eqn{j}-th column is the vector
## of the values of the \eqn{j}-th cardinal basis function on
## the vector \code{xout}.
##
## @param trace Level of verbosity.
##
## @param intOrder Order of the one-dimensional interpolations. Must
## be a permutation of \code{1:d} where \code{d} is the spatial
## dimension.  NOT IMPLEMENTED YET. This argument is similar to the
## argument of the \code{\link{aperm}} method.
## 
## @param ... Further arguments to be passed to \code{interpCB}. NOT
## IMPLEMENTED YET.
##
## @return A single interpolated value if \code{Xout} is a vector or
## a row matrix.  If \code{Xout} is a matrix with several rows, the
## result is a vector of interpolated values, in the order of the
## rows of \code{Xout}.
## 
## @author Yves Deville
##
## @examples
## ## Natural spline for use through Cardinal Basis in 'gridIntCB'
## myInterpCB <-  function(x, xout) cardinalBasis_natSpline(x = x, xout = xout)$CB
## ## Natural spline for use through Cardinal Basis
## myInterp <- function(x, y, xout) {
##      spline(x = x, y = y, n = 3 * length(x), method = "natural", xout = xout)$y
## }
##
## ## generate Grid and function
## set.seed(2468)
## d <- 5
## nLev <- 4L + rpois(d, lambda = 4)
## a <- runif(d)
## myFun2 <- function(x) exp(-crossprod(a, x^2))
## myGD2 <- Grid(nlevels = nLev)
## Y2 <- apply_Grid(myGD2, myFun2)
## n <- 10
## Xout3 <- matrix(runif(n * d), ncol = d)
## t1 <- system.time(GI1 <- gridInt(X = myGD2, Y = Y2, Xout = Xout3, interpFun = myInterp))
## t2 <- system.time(GI2 <- gridInt(X = myGD2, Y = Y2, Xout = Xout3, interpFun = myInterp,
##                                  useC = FALSE))
## t3 <- system.time(GI3 <- gridIntCB(X = myGD2, Y = Y2, Xout = Xout3, interpCB = myInterpCB))
## df <- data.frame(true = apply(Xout3, 1, myFun2),
##                  gridInt_C = GI1, gridInt_R = GI2, gridIntCB = GI3)
## head(df)
## rbind(gridInt_C = t1, gridInt_C = t2, gridIntCB = t3)
##
gridIntCB <- function(X, Y, Xout,
                      interpCB = function(x, xout){ cardinalBasis_lagrange(x = x, xout = xout)$CB },
                      intOrder = NULL,
                      trace = 1L,
                      ...) {

  ##===========================================================================
  ## Coerce 'X' 
  ##===========================================================================
  if (!is(X, "Grid")) {
      if (!is(X, "matrix") && !is(X, "data.frame")) {
          stop("'X' must be a 'Grid' object or a two dimensional",
               " object matrix or data.frame")
      }
      Xco <- TRUE
      X <- as.Grid(X)
      Y <- Y[X@index]
  }
  
  xLevels <- smint::levels(X)
  cxLevels <- lapply(xLevels, as.character)
  d <- dim(X)
  nLevels <- as.integer(smint::nlevels(X))
  ## print(xLevels)
  ## print(nLevels)
  rx <- sapply(xLevels, range)
  nNodes <- prod(nLevels)
  
  ##===========================================================================
  ## check that 'Xout' is correct:
  ##    o a vector of length 'd',
  ##    o a matrix with d columns.
  ## XXX Note that in the future, interpolating on a new (finer) Grid
  ## could be implemented.
  ##===========================================================================
  if (is.matrix(Xout)) {
    if (ncol(Xout) != d) stop("matrix 'X' with incorrect number of columns")
    nOut <- nrow(Xout)
  } else {
    if (length(Xout) != d) stop("when 'X' is not a matrix, it must be a",
                " vector of length 'd', the interpolating dim")
    Xout <- matrix(Xout, nrow = 1)
    nOut <- 1L
  }
  if (trace) {
    cat(sprintf("number of nodes :               %d\n", nNodes))
    cat(sprintf("number of interpolated values:  %d\n", nOut))
  }
  ## check that each 'Xout' values is in the interpolation range 
  rXout <- apply(Xout, 2, range)
  ind <- (rXout[1L, ] < rx[1L, ])
  if (any(ind)) {
    stop("'Xout' values too small for cols ", (1:d)[ind])
  }
  ind <- (rXout[2L, ] > rx[2L, ])
  if (any(ind)) {
    stop("'Xout' values too large for cols ", (1:d)[ind])
  }
  
  ##===========================================================================
  ## check that 'Y' is correct. Accepted objects are
  ##    o a vector with length 'nNodes'.
  ##    o a matrix with 'nNodes' rows and 1 column <-> vector
  ##    o an array 
  ## 
  ##    XXX: list, function ??? NOT IMPLEMENTED YET
  ##
  ##===========================================================================

  if (length(Y) != nNodes) stop("length(Y) must be equal to the number of nodes")
  if (is.array(Y)) dim(Y) <- NULL
  
  ##===========================================================================
  ## Check 'interpCB'. Formals must include 'x' and 'xout'. They will be
  ## reordered if necessary to match this order.
  ##===========================================================================
  if (is.character(interpCB)) interpCB <- get(interpCB, mode = "function")
  fm <- formals(interpCB) 
  if (!all(names(fm)[1L:2L]  == c("x", "xout"))) {
    stop("the first 2 formals o the function 'interpCB' must",
         " be \"x\" and \"xout\"") 
  }
  fmNew <- fm
  fmNew[[1L]] <- fmNew[[2L]] <- NULL
  fmNew <- c(list("x" = NULL, "xout" = NULL), fmNew) 
  formals(interpCB) <- fmNew
  
  ##===========================================================================
  ## Array initialisation: first dim is for the response and dim 2, 3, ...
  ## d + 1 correspond to the spatial dimensions 1, 2, ..., d.
  ##===========================================================================

  Yout <- array(Y, dim = nLevels, dimnames = cxLevels)
  nStar <- as.integer(c(1L, cumprod(nLevels)))

  ##===========================================================================
  ## Now pass an array with d dimensions and dimensions c(nOut,
  ## nLevels[1:(d-1)]). Whe now that d >= 2
  ## ===========================================================================
      
  rho <- new.env()
  environment(interpCB) <- rho
  cat("using '.Call'\n")
  ## print(class(xLevels))

  print(interpCB)
  
  res <- .Call("grid_int_CB",
               Yout,
               nLevels,
               nStar,
               xLevels,
               Xout,
               interpCB,
               rho,
               PACKAGE = "smint")
  
  res[1L:nOut]
  
}
