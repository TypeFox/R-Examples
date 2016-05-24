##*****************************************************************************
## Grid interpolation in arbitrary dimension.
##
## The grid interpolation is performed by looping over
## dimensions. For each dimension, a one-dimensional interpolation is
## carried out, leading to a collection of interpolation problems
## each with a dimension reduced by one.
## 
## @title Grid interpolation in arbitrary dimension
##
## @param X An object that can be coerced into \code{Grid}. This
## can be a data frame or a matrix in \emph{Scattered Data style}:
## the number of columns is then equal to the spatial dimension
## \eqn{d}, and the number of rows is then equal to the number of
## nodes \eqn{n}. But it can also be a \code{Grid} object
## previously created.  A data frame or matrix \code{X} will first be
## coerced into \code{Grid} by using the the S3 method
## \code{\link{as.Grid}}.
## 
## @param Y Response to be interpolated. It must be a vector with
## length \eqn{n} equal to the number of nodes. When \code{X}
## is a matrix or a data frame, the elements in \code{Y}
## match rows of \code{X}, so \code{Y[i]} is the value of the
## interpolated function at \code{X[i, ]}. If instead \code{X}
## is an object with call \code{"Grid"} the order elements
## of \code{Y} must be given in the order associated with the
## object \code{X}.
##
## @param Xout Interpolation locations. Can be a vector or a
## matrix. In the first case, the length of the vector must be equal
## to the spatial dimension \eqn{d} as given by the number of columns
## of \code{X} if this is a matrix, or by \code{dim(X)} if \code{X}
## is a \code{Grid} object.  In the second case, each row will be
## considered as a response to be interpolated, and the number of
## columns of \code{Xout} must be equal to the spatial dimension.
##
## @param interpFun The function to interpolate. This function must
## have as its first 3 formals 'x', 'y', and 'xout', as does
## \code{\link{approx}}.  It must also return the vector of
## interpolated values as DOES NOT \code{approx}. In most cases, a
## simple wrapper will be enough to use an existing method of
## interpolation, see \bold{Examples}.
##
## @param intOrder Order of the one-dimensional interpolations. Must
## be a permutation of \code{1:d} where \code{d} is the spatial
## dimension.  By default, the interpolation order is \eqn{d},
## \eqn{d-1}, \dots, \eqn{1}, corresponding to \code{intOrder =
## d:1L}. Note that for the first element of \code{intOrder}, a
## vector of length \code{nrow(Xout)} is passed as \code{xout} formal
## to \code{interpFun}, while the subsequent \code{interFun} calls
## use a \code{xout} of length \code{1}. So the choice of the first
## element can have an impact on the computation time.
##
## @param useC Logical. If \code{TRUE} the computation is done using
## a C program via the \code{.Call}. Otherwise, the computation is
## done entirely in R by repeated use of the \code{apply} function.
## Normally \code{useC = TRUE} should be faster, but it is not always
## the case.
##
## @param trace Level of verbosity.
##
## @param ... Further arguments to be passed to \code{interpFun}.
##
## @return A single interpolated value when \code{Xout} is either a
## vector or a row matrix.  If \code{Xout} is a matrix with several
## rows, the result is a vector of interpolated values, in the order
## of the rows of \code{Xout}.
## 
## @author Yves Deville
##
## @note A future multivariate version to come will allow the
## simultaneous interpolation of several \emph{ responses}. Most
## probably, this possibility will be used by using a different rule
## for the object \code{Y} (matrix or list).
##
## When \code{X} is a \code{Grid} object and \code{Y} is a
## vector, \emph{the user must ensure that the order of nodes is the
## same for inputs and output}. If so, the order of the dimensions in
## \code{X} can be changed using the generalised transposition
## \code{aperm}. This will change the order of the univariate
## interpolations with possible effects on the computation time
## and on the result when the univariate interpolation method is
## not linear (w.r.t. the response).
##
## %% Recall that default rule for flat
## %% format is that \emph{smaller indices vary faster}. Thus, with
## %% inputs \eqn{X_1}{X1}, \eqn{X_2}{X2}, a flat format representation
## %% will contain a first block of lines with \eqn{X_2}{X2} equal to its
## %% first level, while \eqn{X_1}{X1} runs through all its levels.
##
## When \code{X} is a data frame or a matrix and \code{Y} is a
## vector, the values in \code{Y} are matched to the rows of \code{X}
## is the same as the order. This situation typically arises when
## \code{X} and \code{Y} are columns extracted from a same data
## frame.
## 
## @examples
## set.seed(12345)
## ##========================================================================
## ## Select Interpolation Function. This function must have its first 3
## ## formals 'x', 'y', and 'xout', as does 'approx'. It must also return
## ## the vector of interpolated values as DOES NOT 'approx'. So a wrapper
## ## must be xritten.
## ##=======================================================================
## myInterp <- function(x, y, xout) approx(x = x, y = y, xout = xout)$y
##
## ##=======================================================================
## ## ONE interpolation, d = 2. 'Xout' is a vector.
## ##=======================================================================
## myFun1 <- function(x) exp(-x[1]^2 - 3 * x[2]^2)
## myGD1 <- Grid(nlevels = c("X" = 8, "Y" = 12))
## Y1 <- apply_Grid(myGD1, myFun1)
## Xout1 <- runif(2)
## GI1 <- gridInt(X = myGD1,  Y = Y1, Xout = Xout1, interpFun = myInterp)
## c(true = myFun1(Xout1), interp = GI1)
##
## ##=======================================================================
## ## ONE interpolation, d = 7. 'Xout' is a vector.
## ##=======================================================================
## d <- 7; a <- runif(d); myFun2 <- function(x) exp(-crossprod(a, x^2))
## myGD2 <- Grid(nlevels = rep(4L, time = d))
## Y2 <- apply_Grid(myGD2, myFun2)
## Xout2 <- runif(d)
## GI2 <- gridInt(X = myGD2,  Y = Y2, Xout = Xout2, interpFun = myInterp)
## c(true = myFun2(Xout2), interp = GI2)
##
## ##=======================================================================
## ## n interpolations, d = 7. 'Xout' is a matrix. Same grid data and
## ## response as before
## ##=======================================================================
## n <- 30
## Xout3 <- matrix(runif(n * d), ncol = d)
## GI3 <- gridInt(X = myGD2,  Y = Y2, Xout = Xout3, interpFun = myInterp)
## cbind(true = apply(Xout3, 1, myFun2), interp = GI3)
## 
## ##======================================================================
## ## n interpolation, d = 5. 'Xout' is a matrix. Test the effect of the
## ## order of interpolation.
## ##=======================================================================
## d <- 5; a <- runif(d); myFun4 <- function(x) exp(-crossprod(a, x^2))
## myGD4 <- Grid(nlevels = c(3, 4, 5, 2, 6))
## Y4 <- apply_Grid(myGD4, myFun4)
## n <- 100
## Xout4 <- matrix(runif(n * d), ncol = d)
## t4a <- system.time(GI4a <- gridInt(X = myGD4,  Y = Y4, Xout = Xout4,
##                                    interpFun = myInterp))
## t4b <- system.time(GI4b <- gridInt(X = myGD4,  Y = Y4, Xout = Xout4,
##                                    interpFun = myInterp,
##                                    intOrder = 1L:5L))
## cbind(true = apply(Xout4, 1, myFun4), inta = GI4a, intb = GI4b)

gridInt <- function(X, Y, Xout,
                    interpFun = function(x, y, xout) approx(x = x, y = y, xout = xout)$y,
                    intOrder = NULL,
                    useC = TRUE,
                    trace = 1L,
                    ...) {
   
    ##=========================================================================
    ## Coerce 'X' 
    ##=========================================================================
    if (!is(X, "Grid")) {
        if (!is(X, "matrix") && !is(X, "data.frame")) {
            stop("'X' must be a 'Grid' object or a two dimensional",
                 " object matrix or data.frame")
        }
        Xco <- TRUE
        X <- as.Grid(X)
        Y <- Y[X@index]
    }
    d <- dim(X)
    
    ##========================================================================
    ## check that 'Xout' is correct:
    ##    o a vector of length 'd',
    ##    o a matrix with d columns.
    ## XXX Note that in the future, interpolating on a new (finer) Grid
    ## could be implemented.
    ##========================================================================
    if (is.matrix(Xout)) {
        if (ncol(Xout) != d) stop("matrix 'Xout' with incorrect number of columns")
        nOut <- nrow(Xout)
    } else {
        if (length(Xout) != d) stop("when 'Xout' is not a matrix, it must be a",
                      " vector of length 'd', the interpolating dim")
        Xout <- matrix(Xout, nrow = 1)
        nOut <- 1L
    }
    nNodes <- prod(nlevels(X))

    ##========================================================================
    ## check that 'Y' is correct. Accepted objects are
    ##    o a vector with length 'nNodes'.
    ##    o a matrix with 'nNodes' rows and 1 column <-> vector
    ##    o an array 
    ## 
    ##    XXX: list, function ??? NOT IMPLEMENTED YET
    ##
    ##=======================================================================
    if (length(Y) != nNodes) stop("length(Y) must be equal to the number of ",
                  "nodes")
    if (is.array(Y)) dim(Y) <- NULL
    
    ##=========================================================================
    ## general transpose of 'X' if necessary
    ##=========================================================================
    if (is.null(intOrder)) {
        intOrder <- rev(1L:d)
    } else {
        nx <- nlevels(X)
        X <- aperm(X, perm = rev(intOrder))
        Xout <- Xout[ , rev(intOrder), drop = FALSE]
        Y <- array(Y, dim = nx)
        Y <- aperm(Y, perm = rev(intOrder))
    }
    
    xLevels <- levels(X)
    cxLevels <- lapply(xLevels, as.character)
    
    nx <- nlevels(X)
    rx <- sapply(xLevels, range)
     
    if (trace) {
        cat(sprintf("o Number of nodes :               %d\n", nNodes))
        cat(sprintf("o Number of interpolated values:  %d\n", nOut))
        cat("o Interpolation locations")
        if (nOut > 4L)  cat(sprintf(" (first out of %d)\n", nOut))
        else cat("\n")
        print(head(Xout, n = 4L))
    }

    ## check that each 'Xout' values is in the interpolation range 
    rXout <- apply(Xout, 2, range)

    ind <- (rXout[1L, ] < rx[1L, ])
    if (any(ind)) stop("'Xout' values too small for cols ", (1:d)[ind])

    ind <- (rXout[2L, ] > rx[2L, ])
    if (any(ind)) stop("'Xout' values too large for cols ", (1:d)[ind])
    
    ##=======================================================================
    ## Check 'interpfun'. Formal 'y' must now come in first position in
    ## order to use the standard 'apply' facility.
    ##=======================================================================
    if (is.character(interpFun)) interpFun <- get(interpFun, mode = "function")
    fm <- formals(interpFun) 
    if (!all(names(fm)[1L:3L]  == c("x", "y", "xout"))) {
        stop("the first 3 formals of the function 'interpFun' must",
             " be \"x\", \"y\", and \"xout\"") 
    }
    fmNew <- fm
    fmNew[[1L]] <- fmNew[[2L]] <- fmNew[[3L]] <- NULL
    fmNew <- c(list("y" = NULL, "x" = NULL, "xout" = NULL), fmNew) 
    formals(interpFun) <- fmNew
    
    ##========================================================================
    ## Get rid of the one-dimensional case!
    ##========================================================================
    if (d == 1L) {
        Yout <- interpFun(Y, x = xLevels[[1]], xout = Xout[ , 1L])
        return(Yout)
    }
    
    ##========================================================================
    ## Use the .Call interface
    ##========================================================================
    if (useC) {
        
        ## Array initialisation: first dim is for the response and dim 2, 3, ...
        ## d + 1 correspond to the spatial dimensions 1, 2, ..., d.
        Yout <- array(Y, dim = nx, dimnames = cxLevels)
       
        ## first interpolation: after it, the first dimension of the array
        ## has 'nOut' values and the last dimension is lost. 
        Yout <- apply(Yout, MARGIN = 1L:(d - 1L), FUN = interpFun,
                      x = xLevels[[d]], xout = Xout[ , d])
        nStar <- as.integer(c(1L, cumprod(nx)))
        
        ## Now pass an array with d dimensions and dimensions c(nOut,
        ## nLevels[1:(d-1)]). Whe know that d >= 2
        rho <- new.env()
        environment(interpFun) <- rho
        if (trace) cat("o Using '.Call'\n")
        
        res <- .Call("grid_int",
                     Yout,
                     nx,
                     nStar,
                     xLevels,
                     Xout,
                     interpFun,
                     rho,
                     PACKAGE = "smint")
        
        return(res[1L:nOut])

    }
    
    ##========================================================================
    ## Array initialisation
    ##========================================================================
    if (any(nx == 1)) 
        interpFun1 <- function(x, y, xout) rep(y[1], length.out = length(xout))
    
    Y <- array(Y, dim = nx, dimnames = cxLevels)

    ##========================================================================
    ## When nOut == 1, the code is simpler
    ##========================================================================
    if (nOut == 1L) {
        ## add one dimension
        Y <- array(Y, dim = c(1, nx), dimnames = c("out", dimnames(Y)))
        ## undidimensional interpolation in dimensions d + 1, d, ..., 2
        d1 <- d + 1L
        dj <- d1
        dimj <- dim(Y)
        
        for (j in d1:2L) {
            dj <- dj - 1L
            if (nx[j - 1L] > 1L) intFun <- interpFun
            else intFun <- interpFun1
            Y <- apply(Y, MARGIN = 1L:dj, FUN = intFun,
                       x = xLevels[[j - 1L]],
                       xout = Xout[ , j - 1L])
            dimj[j] <- 1L
            dim(Y) <- dimj
            if (trace) {
                cat(sprintf("   interp. %2d,  var. %8s, nx = %3d, nxout = %3d\n",
                            d1 - j + 1, names(xLevels)[j - 1L],
                            length(xLevels[[j - 1L]]),
                            nrow(Xout)))
                cat("levels = ", xLevels[[j - 1L]], "xout = ", Xout[ , j - 1L], "\n")
            }
        }
        
        return(Y)
    }
    
    ##========================================================================
    ## first interpolation: after it, the first dimension of the array
    ## matches the rows of 'Xout': it has 'nOut' values and the last
    ## dimension is lost.
    ## =======================================================================
    if (nx[d] > 1L) intFun <- interpFun
    else intFun <- interpFun1
    
    Y <- apply(Y, MARGIN = 1L:(d - 1L), FUN = intFun, x = xLevels[[d]],
               xout = Xout[ , d])

    if (trace) {
        cat(sprintf("   interp. %2d,  var. %8s, nx = %3d, nxout = %3d\n",
                    1, names(xLevels)[d],
                    length(xLevels[[d]]), nrow(Xout)))
    }

    ##========================================================================
    ## Intermediate interpolations: after each of them, the first
    ## dimension of the array has 'nOut' values and the last dimension
    ## is lost.
    ## =======================================================================
    if (d > 2L) {
        for (j in (d - 1L):2L) {
            Yout2 <- array(NA, dim = c(nOut, nx[1L:(j - 1L)]))
            if (nx[j] > 1L) intFun <- interpFun
            else intFun <- interpFun1
            for (k in 1L:nOut) {
                Youtk <- Y[slice.index(Y, MARGIN = 1L) == k] ## 'j' dimensions
                dim(Youtk) <- nx[1L:j]
                Youtk <- apply(Youtk, MARGIN = 1L:(j - 1L), FUN = intFun,
                               x = xLevels[[j]], xout = Xout[k, j])
                Yout2[slice.index(Yout2, MARGIN = 1L) == k] <- Youtk
            }
            Y <- Yout2
            if (trace) {
                cat(sprintf("   interp. %2d,  var. %8s, nx = %3d, nxout = %3d\n",
                            d - j + 1L, names(xLevels)[j],
                            length(xLevels[[j]]), nOut))
            }
        }
        
    } 
    ## d == 2: back to a one-dimensional interpolation
    j <- 1L
    Yout2 <- rep(NA, nOut)
    if (nx[1] > 1L) intFun <- interpFun
    else intFun <- interpFun1
    for (k in 1L:nOut) {
        Yout2[k] <- intFun(Y[k, ], x = xLevels[[1L]], xout = Xout[k, 1])
    }
    if (trace) {
        cat(sprintf("   interp. %2d,  var. %8s, nx = %3d, nxout = %3d\n",
                    d - j + 1L, names(xLevels)[j],
                    length(xLevels[[j]]), nOut))
            }
    return(Yout2)
    
}
