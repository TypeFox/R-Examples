##' Grid linear interpolation in arbitrary dimension through repeated
##' one-dimensional interpolations.
##'
##' The one-dimensional interpolation is performed using an
##' interpolation function given in \code{interpFun1d} or a Cardinal
##' Basis function given in \code{cardinalBasis1d}.  In the second
##' case, it is required here that the one-dimensional interpolation
##' method is \emph{linear} w.r.t. the vector of interpolated values.
##' For each dimension, a one-dimensional interpolation is carried
##' out, leading to a collection of interpolation problems each with a
##' dimension reduced by one.
##'
##' When a cardinal basis is used, the \emph{same basis} is used for
##' all interpolations w.r.t. the same variable, so the number of call
##' to the \code{cardinalBasis1d} function is equal to the interpolation
##' dimension, making this method \emph{very fast} compared to the
##' general grid interpolation, see the \bold{Examples} section.
##'
##' By default, i.e. when both \code{interpFun1d} and \code{cardinalBasis1d}
##' are \code{NULL}, a Lagrange cardinal basis interpolation is used. So
##' the 1d interpolation is the broken line interpolation.
##' 
##' @title Grid interpolation in arbitrary dimension through Cardinal
##' Basis
##'
##' @aliases interp_Grid
##' 
##' @usage
##' interp_Grid(X, Y, Xout,
##'             interpFun1d = NULL,
##'             cardinalBasis1d = NULL,
##'             intOrder = NULL,
##'             useC = TRUE, trace = 1L,
##'             ...)
##' 
##' @param X An object that can be coerced into \code{Grid}. This
##' can be a data.frame or a matrix in \emph{Scattered Data style} in
##' which case the column number is equal to the spatial dimension
##' \eqn{d}, and the row number is then equal to the number of nodes
##' \eqn{n}. But it can also be a \code{Grid} object previously
##' created.  A data frame or matrix \code{X} will first be coerced
##' into \code{Grid} by using the the S3 method
##' \code{\link{as.Grid}}.
##' 
##' @param Y Response to be interpolated. It must be a vector of
##' length \eqn{n} equal to the number of nodes. When \code{X} has
##' class \code{Grid}, the order of the elements in \code{Y} must
##' conform to the order of the nodes as given in \code{X}.
##'
##' @param Xout Interpolation locations. Can be a vector or a
##' matrix. In the first case, the length of the vector must be equal
##' to the spatial dimension \eqn{d} as given by \code{xLev} or
##' \code{X}.  In the second case, each row will be considered as a
##' response to be interpolated, and the number of columns of
##' \code{Xout} must be equal to the spatial dimension.
##'
##' @param interpFun1d The function to interpolate in 1d. This
##' function must have as its first 3 formals 'x', 'y', and 'xout', as
##' does \code{\link{approx}}.  It must also return the vector of
##' interpolated values as DOES NOT \code{approx}. In most cases, a
##' simple wrapper will be enough to use an existing method of
##' interpolation, see \bold{Examples}.
##' 
##' @param cardinalBasis1d Function evaluating an interpolation
##' Cardinal Basis. This function must have as its first 2 formals
##' 'x', and 'xout'.  It must return a matrix with \code{length(x)}
##' columns and \code{length(xout)} rows. The \eqn{j}-th column is the
##' vector of the values of the \eqn{j}-th cardinal basis function on
##' the vector \code{xout}. \emph{This formal argument can be used
##' only when} \code{interpFun1d} \emph{is not used}.
##'
##' @param trace Level of verbosity.
##'
##' @param intOrder Order of the one-dimensional interpolations.  Used
##' only when \code{interpFun1d} is given, and not when
##' \code{cardinalBasis1d} is given. Must be a permutation of
##' \code{1:d} where \code{d} is the spatial dimension. This argument
##' is similar to the argument of the \code{\link{aperm}} method.  By
##' default, the interpolation order is \eqn{d}, \eqn{d-1}, \dots,
##' \eqn{1}, corresponding to \code{intOrder = d:1L}. Note that for
##' the first element of \code{intOrder}, a vector of length
##' \code{nrow(Xout)} is passed as \code{xout} formal to
##' \code{interpFun1d}, while the subsequent calls to
##' \code{interpFun1} use a \code{xout} of length \code{1}. So the
##' choice of the first element can have an impact on the computation
##' time.
##' 
##' @param useC Logical. Used only when \code{interpFun1d} is given,
##' and not when \code{cardinalBasis1d} is given. If \code{TRUE}, the
##' computation is done using ' a C program via the
##' \code{.Call}. Otherwise, the computation is ' done entirely in R
##' by repeated use of the \code{apply} function.  ' Normally
##' \code{useC = TRUE} should be faster, but it is not always ' the
##' case.
##'
##' @param ... Further arguments to be passed to \code{interpCB}. NOT
##' IMPLEMENTED YET.
##'
##' @return A single interpolated value if \code{Xout} is a vector or
##' a row matrix.  If \code{Xout} is a matrix with several rows, the
##' result is a vector of interpolated values, in the order of the
##' rows of \code{Xout}.
##' 
##' @author Yves Deville
##'
##' @examples
##' set.seed(12345)
##' ##========================================================================
##' ## Select Interpolation Function. This function must have its first 3
##' ## formals 'x', 'y', and 'xout', as does 'approx'. It must also return
##' ## the vector of interpolated values as DOES NOT 'approx'. So a wrapper
##' ## must be xritten.
##' ##=======================================================================
##' myInterpFun <- function(x, y, xout) approx(x = x, y = y, xout = xout)$y
##'
##' ##=======================================================================
##' ## Example 1
##' ## ONE interpolation, d = 2. 'Xout' is a vector.
##' ##=======================================================================
##' myFun1 <- function(x) exp(-x[1]^2 - 3 * x[2]^2)
##' myGD1 <- Grid(nlevels = c("X" = 8, "Y" = 12))
##' Y1 <- apply_Grid(myGD1, myFun1)
##' Xout1 <- runif(2)
##' GI1 <- interp_Grid(X = myGD1,  Y = Y1, Xout = Xout1,
##'                    interpFun1d = myInterpFun)
##' c(true = myFun1(Xout1), interp = GI1)
##'
##' ##=======================================================================
##' ## Example 2
##' ## ONE interpolation, d = 7. 'Xout' is a vector.
##' ##=======================================================================
##' d <- 7; a <- runif(d); myFun2 <- function(x) exp(-crossprod(a, x^2))
##' myGD2 <- Grid(nlevels = rep(4L, time = d))
##' Y2 <- apply_Grid(myGD2, myFun2)
##' Xout2 <- runif(d)
##' GI2 <- interp_Grid(X = myGD2,  Y = Y2, Xout = Xout2,
##'                    interpFun1d = myInterpFun)
##' c(true = myFun2(Xout2), interp = GI2)
##'
##' ##=======================================================================
##' ## Example 3
##' ## 'n' interpolations, d = 7. 'Xout' is a matrix. Same grid data and
##' ## response as before
##' ##=======================================================================
##' n <- 30
##' Xout3 <- matrix(runif(n * d), ncol = d)
##' GI3 <- interp_Grid(X = myGD2,  Y = Y2, Xout = Xout3,
##'                    interpFun1d = myInterpFun)
##' cbind(true = apply(Xout3, 1, myFun2), interp = GI3)
##' 
##' ##======================================================================
##' ## Example 4
##' ## 'n' interpolation, d = 5. 'Xout' is a matrix. Test the effect of the
##' ## order of interpolation.
##' ##=======================================================================
##' d <- 5; a <- runif(d); myFun4 <- function(x) exp(-crossprod(a, x^2))
##' myGD4 <- Grid(nlevels = c(3, 4, 5, 2, 6))
##' Y4 <- apply_Grid(myGD4, myFun4)
##' n <- 100
##' Xout4 <- matrix(runif(n * d), ncol = d)
##' t4a <- system.time(GI4a <- interp_Grid(X = myGD4,  Y = Y4, Xout = Xout4,
##'                                        interpFun1d = myInterpFun))
##' t4b <- system.time(GI4b <- interp_Grid(X = myGD4,  Y = Y4, Xout = Xout4,
##'                                        interpFun1d = myInterpFun,
##'                                        intOrder = 1L:5L))
##' cbind(true = apply(Xout4, 1, myFun4), inta = GI4a, intb = GI4b)
##'
##' ##======================================================================
##' ## Example 5 
##' ## 'n' interpolation, d = 5 using a GIVEN CARDINAL BASIS.
##' ## 'Xout' is a matrix. 
##' ##=======================================================================
##'
##' ## Natural spline for use through Cardinal Basis in 'gridIntCB'
##' myInterpCB5 <-  function(x, xout) cardinalBasis_natSpline(x = x, xout = xout)$CB
##' 
##' ## Natural spline for use through an interpolation function
##' myInterpFun5 <- function(x, y, xout) {
##'      spline(x = x, y = y, n = 3 * length(x), method = "natural", xout = xout)$y
##' }
##' 
##' n <- 10
##' ## Since we use a natural spline, there must be at least 3 levels for
##' ## each dimension.
##' myGD5 <- Grid(nlevels = c(3, 4, 5, 4, 6))
##' Y5 <- apply_Grid(myGD5, myFun4)
##' Xout5 <- matrix(runif(n * d), ncol = d)
##' 
##' t5a <- system.time(GI5a <- interp_Grid(X = myGD5, Y = Y5, Xout = Xout5,
##'                                        interpFun1d = myInterpFun5))
##' t5b <- system.time(GI5b <- interp_Grid(X = myGD5, Y = Y5, Xout = Xout5,
##'                                        interpFun1d = myInterpFun5,
##'                                        useC = FALSE))
##' t5c <- system.time(GI5c <- interp_Grid(X = myGD5, Y = Y5, Xout = Xout5,
##'                                        cardinalBasis1d = myInterpCB5))
##' df <- data.frame(true = apply(Xout5, 1, myFun4),
##'                  gridInt_C = GI5a, gridInt_R = GI5b, gridIntCB = GI5c)
##' head(df)
##' rbind(gridInt_C = t5a, gridInt_R = t5b, gridIntCB = t5c)
##' 
interp_Grid <- function(X, Y, Xout,
                        interpFun1d = NULL,
                        cardinalBasis1d = NULL,
                        intOrder = NULL,
                        useC = TRUE,
                        trace = 1L,
                        ...) {

    if (is.null(cardinalBasis1d)) {
        if (is.null(interpFun1d)) {
            type <- "CB"
            cardinalBasis1d <- function(x, xout){
                cardinalBasis_lagrange(x = x, xout = xout)$CB
            }
        } else {
            type <- "fun"
        }
    } else {
        if (!is.null(interpFun1d)) {
            stop("only one of the formals 'interpFun1d' and",
                 " 'cardinalBasis1s' can be given")
        }
        type <- "CB"
    }

    if (type == "CB") {
        res <- gridIntCB(X, Y, Xout,
                         interpCB = cardinalBasis1d,
                         intOrder = intOrder,
                         trace = 1L,
                         ...) 
    } else {
        res <- gridInt(X, Y, Xout,
                       interpFun = interpFun1d,
                       intOrder = intOrder,
                       useC = useC,
                       trace = 1L,
                       ...) 
    }
    
    res
    
}
