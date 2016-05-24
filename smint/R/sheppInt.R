
##*****************************************************************************
##' Shepard's (modified) quadratic interpolation method. 
##'
##' Shepard's modified interpolation method works for scattered data
##' in arbitrary dimension. It relies on a \emph{local} polynomial fit
##' and the quadratic version uses a polynomial of degree 2 (a quadratic
##' form) as local approximation of the function.
##' 
##' @title Shepard's (modified) quadratic interpolation method 
##'
##' @aliases qsheppInt qsheppInt2d qsheppInt3d
##'
##' @usage
##' 
##' qsheppInt(X, y, XNew = NULL, nQ = 15L, nW = 33L, nR,
##'           checkX = TRUE) 
##' qsheppInt2d(X, y, XNew = NULL, nQ = 13L, nW = 19L, nR,
##'             deriv = 0L, checkX = TRUE)
##' qsheppInt3d(X, y, XNew = NULL, nQ = 17L, nW = 32L, nR,
##'             deriv = 0L, checkX = TRUE) 
##' 
##' @param X Design matrix. A matrix with \eqn{n} rows and \eqn{d}
##' columns where \eqn{n} is the number of nodes and \eqn{d} is the
##' dimension.
##'
##' @param y Numeric vector of length \eqn{n} containing the function
##' values at design points.
##' 
##' @param XNew Numeric matrix of new design points where
##' interpolation is computed. It must have \eqn{d} columns.
##'
##' @param nQ Number of points used in the local polynomial fit. This
##' is the parameter \code{NQ} of the Fortran routine \code{qshepmd}
##' and the parameter \eqn{N_p} of the referenced article.
##'
##' @param nW Number of nodes within (and defining) the radii of
##' influence, which enter into the weights. This is parameter
##' \code{NW} of the Fortran routine \code{qshepmd} and \eqn{N_w} of
##' the reference article.
##'
##' @param nR Number of divisions in each dimension for the cell grid
##' defined in the subroutine \code{storem}.  A hyperbox containing
##' the nodes is partitioned into cells in order to increase search
##' efficiency.  The recommended value \eqn{(n/3)^(1/d)} is used as
##' default, where \eqn{n} is the number of nodes and \eqn{d} is the
##' dimension.
##'
##' @param deriv Logical or integer value in \code{c(0, 1)}. When the
##' (coerced) integer value is \code{1}, the derivatives are computed.
##' This is is only possible for the dimensions \code{2} and \code{3},
##' and only when the specific interpolation functions are used.
##' The result is then a matrix with one column for the function and
##' one column by derivative.
##' 
##' @param checkX If \code{TRUE}, the dimensions and colnames of
##' \code{X} and \code{XNew} will be checked.
##'
##' @return A list with several objects related to the computation
##' method. The vector of interpolated value is in the list element
##' named \code{yNew}.
##' 
##' @note This function is an R interface to the \code{qshepmd}
##' routine in the \pkg{SHEPPACK} Fortran package available on netlib
##' \url{http://www.netlib.org} as algorithm 905A.
##'
##' The \code{qshepInt} function is an interface for the
##' \code{QSHEPMD} Fortran routine, while \code{qshepInt2d} and
##' \code{qshepInt3d} are interfaces to the \code{QSHEPM2D} and
##' \code{QSHEPM3D} Fortran routines. The general interpolation of
##' \code{qshepInt} can be used also for the dimensions \eqn{2} and
##' \eqn{3}. However, this function does not allow the computation of
##' the derivatives as \code{qshepInt2d} and \code{qshepInt3d} do.
##' 
##' @author Fortran code by William I. Thacker; Jingwei
##' Zhang; Layne T. Watson; Jeffrey B. Birch; Manjula A. Iyer; Michael
##' W. Berry. See References below.  The Fortran code is a translation
##' of M.W Berry's C++ code.
##' 
##' Code adaptation and R interface by Yves Deville.
##' 
##' @references W.I. Thacker, J. Zhang, L.T. Watson, J.B. Birch,
##' M.A. Iyer and M.W. Berry (2010). Algorithm 905: SHEPPACK: Modified
##' Shepard Algorithm for Interpolation of Scattered Multivariate Data
##' \emph{ACM Trans. on Math. Software} (TOMS) Vol. 37, n. 3.
##' \href{http://dl.acm.org/citation.cfm?id=1824812}{link}
##' 
##' M.W. Berry and K.S. Minser (1999). Algorithm 798: High-dimensional
##' interpolation using the modified Shepard method. \emph{ACM
##' Trans. Math. Software} (TOMS) Vol. 25, n. 3, pp. 353-366.
##' \href{http://dl.acm.org/citation.cfm?id=326147.326154}{link}
##' 
##' @examples
##' n <- 1500; nNew <- 100; d <- 4
##' fTest <- function(x)((x[1] + 2 * x[2]   + 3 * x[3] + 4 * x[4]) / 12)^2
##' set.seed(12345)
##' X <- matrix(runif(n*d), nrow = n, ncol = d)
##' y <- apply(X, 1, FUN = fTest)
##' XNew <- matrix(runif(nNew * d), nrow = nNew, ncol = d)
##' yNew <- apply(XNew, 1, FUN = fTest)
##' system.time(res <- qsheppInt(X = X, XNew = XNew, y = y, nQ = 40,
##'                              checkX = FALSE))
##' ## check errors
##' max(abs(res$yNew - yNew))
##' 
##' ##=========================================================================
##' ## Use SHEPPACK test functions see Thacker et al. section 7 'PERFORMANCE'
##' ##=========================================================================
##' \dontrun{
##'    set.seed(1234)
##'    d <- 3
##'    k <- 0:4; n0 <- 100 * 2^k; n1 <- 4
##'    GD <- Grid(nlevels = rep(n1, d))
##'    XNew <- as.matrix(GD)
##'    RMSE <- array(NA, dim = c(5, length(k)),
##'                  dimnames = list(fun = 1:5, k = k))
##'    for (iFun in 1:5) {
##'       yNew <- apply(XNew, 1, ShepFuns[[iFun]])
##'       for (iN in 1:length(n0)) {
##'          X <- matrix(runif(n0[iN] * d), ncol = d)
##'          y <- apply(X, 1, ShepFuns[[iFun]])
##'          res <- qsheppInt(X = X, XNew = XNew, y = y, nQ = 40, checkX = FALSE)
##'          RMSE[iFun, iN] <- mean((res$yNew - yNew)^2)
##'       }
##'    }
##'    cols <- c("black", "SteelBlue2", "orangered", "SpringGreen3", "purple")
##'    pchs <-  c(16, 21, 22, 23, 24)
##'    matplot(k, t(RMSE), type = "o", lwd = 2, lty = 1,
##'            col = cols, xaxt = "n", pch = pchs, cex = 1.4,
##'            bg = "white",
##'            main = sprintf("dim = %d SHEPPACK test functions", d),
##'            xlab = "number of nodes", ylab = "RMSE")
##'    axis(side = 1, at = k, labels = n0)
##'    legend("topright", legend = paste("shepFun", 1:5),
##'            col = cols, pch = pchs, lwd = 2, pt.lwd = 2, pt.bg = "white")
##' }
qsheppInt <- function(X, y, XNew = NULL,
                      nQ = 15L, nW = 33L, nR,
                      checkX = TRUE) {

    if (nQ < 1) stop ("found 'nQ < 1'")
    
    ## check matrices
    if (checkX) {
        XNew <- checkX(X = X, XNew = XNew)
        X <- XNew$X
        XNew <- XNew$XNew
    }
    
    d <- ncol(X)
    n <- nrow(X)
    nNew <- nrow(XNew)
    yNew <- rep(0, nNew)
    
    NTT <- (d * (d + 3L)) / 2 + 1

    ## intercept easily detected errors before calling Fortran
    if (nQ < NTT -1) stop("found 'nQ < NTT - 1'")

    if (missing(nR)) nR <- ceiling((n / 3)^(1 / d))

    if (nR < 1) stop("found 'nR < 1'")
    
    NQWMAX <- pmax(nQ, nW)
    LMAX <- pmin(100, n - 1L)
    if (NQWMAX > LMAX) {
        stop("'nQ' and 'nW' must be <= 100 and <= n")
    }

    ## allocate space from R, since it is not possible to do this
    ## in Fortran
    A <- numeric(n * NTT)
    DX <- numeric(d)
    XMIN <- numeric(d)
    RSQ <- numeric(n)
    LCELL <- integer(nR^d)
    LNEXT <- integer(n)
    WS <- numeric(NTT*NTT)
    IW <- numeric(d*5)
    
    ## call Fortran
    res <- .Fortran("SHEPQM",
                    d = as.integer(d),         ## M
                    n = as.integer(n),         ## N
                    X = as.double(t(X)),       ## X
                    y = as.double(y),          ## F
                    nNew = as.integer(nNew),   ## NNEW
                    XNew = as.double(t(XNew)), ## XNEW
                    yNew = as.double(yNew),    ## FNEW
                    nQ = as.integer(nQ),       ## NQ
                    nW = as.integer(nW),       ## NW
                    nR = as.integer(nR),       ## NR
                    rMax = as.double(1.0),     ## RMAX
                    IER = as.integer(0L),      ## IER
                    ## new arguments to avoid allocation (from the stack)
                    NTT = as.integer(NTT),
                    A = as.double(A),
                    DX = as.double(DX),
                    XMIN = as.double(XMIN),
                    RSQ = as.double(RSQ),
                    LCELL = as.integer(LCELL),
                    LNEXT = as.integer(LNEXT),
                    WS = as.double(WS),
                    IW = as.double(IW),
                    PACKAGE = "smint")
       
    ## back to matrix 'X' 
    res$X <- X
    res$XNew <- XNew
    
    ## manage error codes
    if (res$IER == 1L) {
        stop("'IER = 1'. 'M', 'N', 'NQ', 'NW', or 'NR' is out of range")
    }
    if (res$IER == 2L) stop("'IER = 2'. Duplicate nodes were encountered")
    if (res$IER == 3L) {
        stop("'IER = 3'. All the nodes lie in an affine subspace",
             " of dimension M-1")
    }
    if (res$IER == 4L) {
        stop("'IER = 4'. Space could not be allocated for",
             " temporary arrays")
    }
    return(res)
    
}

qsheppInt2d <- function(X, y, XNew = NULL,
                        nQ = 13L, nW = 19L, nR,
                        deriv = 0L,
                        checkX = TRUE) {
    
    if (nQ < 1) stop ("found 'nQ < 1'")

    deriv <- as.integer(deriv)
    
    if (deriv == 0) {
        nc <- 1L
        cNm <- c("y")
        routineName <- "SHEPQ2"
    } else if (deriv == 1) {
        nc <- 3L
        cNm <- c("y", "dX1", "dX2")
        routineName <- "SHEPQ2G"
    } else {
        stop("'deriv' must be equal to 0 or 1")
    }
    
    ## check matrices
    if (checkX) {
        XNew <- checkX(X = X, XNew = XNew)
        X <- XNew$X
        XNew <- XNew$XNew
    }
    d <- ncol(X)
    n <- nrow(X)
    if (d != 2L) stop("this method is only for dimension 2")
    if (length(y) != n) stop("bad length for 'y'")
   
    nNew <- nrow(XNew)
    yNew <- rep(0, nNew * nc)

    if (missing(nR)) nR <- ceiling((n / 3)^(1 / d))
    
    ## intercept easily detected errors before calling Fortran
    if (nQ < 5) stop("found 'nQ < 5'")
    if (nR < 1) stop("found 'nR < 1'")
    if (nW < 1) stop("found 'nW < 1'")
    
    NQWMAX <- pmax(nQ, nW)
    LMAX <- pmin(40, n - 1L)
    if (NQWMAX > LMAX) {
        stop("'nQ' and 'nW' must be <= 40 and <= n ")
    }

    ## these variables do not need initialisation, just
    ## storage
    LCELL <- integer(nR^2)
    LNEXT <- integer(n)
    XMIN <- numeric(1L)
    YMIN <- numeric(1L)
    DX <- numeric(1L)
    DY <- numeric(1L)
    RSQ <- numeric(n)
    A <- numeric(n * 5L)
    
    ## call Fortran
    res <- .Fortran(routineName,
                    n = as.integer(n),         ## N
                    X = as.double(X[ , 1L]),   ## X
                    Y = as.double(X[ , 2L]),   ## Y
                    y = as.double(y),          ## F
                    nNew = as.integer(nNew),   ## NNEW
                    XNew = as.double(XNew[ , 1L]),    ## XNEW
                    YNew = as.double(XNew[ , 2L]),    ## XYEW
                    yNew = as.double(yNew),    ## FNEW
                    nQ = as.integer(nQ),       ## NQ
                    nW = as.integer(nW),       ## NW
                    nR = as.integer(nR),       ## NR
                    LCELL = as.integer(LCELL),
                    LNEXT = as.integer(LNEXT),
                    XMIN = as.double(XMIN),
                    YMIN = as.double(YMIN),
                    DX = as.double(DX),
                    DY = as.double(DY),
                    rMax = as.double(1.0),    ## RMAX
                    RSQ = as.double(RSQ),
                    A = as.double(A),
                    IER = as.integer(0L),      ## IER
                    PACKAGE = "smint")

    res$yNew <- matrix(res$yNew, nrow = nNew, ncol = nc)
    colnames(res$yNew) <- cNm
    
    ## back to matrix 'X' 
    res$X <- X
    res$Y <- NULL
    res$XNew <- XNew
    res$YNew <- NULL
    
    ## manage error codes
    if (res$IER == 1L) {
        stop("'IER = 1'. 'M', 'N', 'NQ', 'NW', or 'NR' is out of range")
    }
    if (res$IER == 2L) stop("'IER = 2'. Duplicate nodes were encountered")
    if (res$IER == 3L) {
        stop("'IER = 3'. All the nodes lie in an affine subspace",
             " of dimension M-1")
    }

    return(res)
    
}

qsheppInt3d <- function(X, y, XNew = NULL,
                        nQ = 17L, nW = 32L, nR,
                        deriv = 0L,
                        checkX = TRUE) {
    
    if (nQ < 1) stop ("found 'nQ < 1'")

    deriv <- as.integer(deriv)
    
    if (deriv == 0) {
        nc <- 1L
        cNm <- c("y")
        routineName <- "SHEPQ3"
    } else if (deriv == 1) {
        nc <- 4L
        cNm <- c("y", "dX1", "dX2", "dX3")
        routineName <- "SHEPQ3G"
    } else {
        stop("'deriv' must be equal to 0 or 1")
    }
    
    ## check matrices
    if (checkX) {
        XNew <- checkX(X = X, XNew = XNew)
        X <- XNew$X
        XNew <- XNew$XNew
    }
    d <- ncol(X)
    n <- nrow(X)
    if (d != 3L) stop("this method is only for dimension 3")
    if (length(y) != n) stop("bad length for 'y'")
    
    nNew <- nrow(XNew)
    yNew <- rep(0, nNew * nc)

    if (missing(nR)) nR <- ceiling((n / 3)^(1 / d))
    
    ## intercept easily detected errors before calling Fortran
    if (nQ < 9) stop("found 'nQ < 9'")
    if (nR < 1) stop("found 'nR < 1'")
    if (nW < 1) stop("found 'nW < 1'")
    
    NQWMAX <- pmax(nQ, nW)
    LMAX <- pmin(40, n - 1L)
    if (NQWMAX > LMAX) {
        stop("'nQ' and 'nW' must be <= 40 and <= n ")
    }

    ## these variables do not need initialisation, just
    ## storage
    LCELL <- integer(nR^3)
    LNEXT <- integer(n)
    XYZMIN <- numeric(3L)
    XYZDEL <- numeric(3L)
    RSQ <- numeric(n)
    A <- numeric(n * 9L)
    
    ## call Fortran
    res <- .Fortran(routineName,
                    n = as.integer(n),              ## N
                    X = as.double(X[ , 1L]),        ## X
                    Y = as.double(X[ , 2L]),        ## Y
                    Z = as.double(X[ , 3L]),        ## Z
                    y = as.double(y),               ## F
                    nNew = as.integer(nNew),        ## NNEW
                    XNew = as.double(XNew[ , 1L]),  ## XNEW
                    YNew = as.double(XNew[ , 2L]),  ## YNEW
                    ZNew = as.double(XNew[ , 3L]),  ## ZNEW
                    yNew = as.double(yNew),         ## FNEW
                    nQ = as.integer(nQ),            ## NQ
                    nW = as.integer(nW),            ## NW
                    nR = as.integer(nR),            ## NR
                    LCELL = as.integer(LCELL),   
                    LNEXT = as.integer(LNEXT),
                    XYZMIN = as.double(XYZMIN),
                    XYZDEL = as.double(XYZDEL),
                    rMax = as.double(1.0),          ## RMAX
                    RSQ = as.double(RSQ),           ## RSQ
                    A = as.double(A),               ## A
                    IER = as.integer(0L),           ## IER
                    PACKAGE = "smint")

    res$yNew <- matrix(res$yNew, nrow = nNew, ncol = nc)
    colnames(res$yNew) <- cNm

    ## back to matrix 'X' 
    res$X <- X
    res$Y <- res$Z <- NULL
    res$XNew <- XNew
    res$YNew <- res$ZNew <- NULL
    
    ## manage error codes
    if (res$IER == 1L) {
        stop("'IER = 1'. 'M', 'N', 'NQ', 'NW', or 'NR' is out of range")
    }
    if (res$IER == 2L) stop("'IER = 2'. Duplicate nodes were encountered")
    if (res$IER == 3L) {
        stop("'IER = 3'. All the nodes lie in an affine subspace",
             " of dimension M-1")
    }

    return(res)
    
}


## lsheppIntmd <- function(X, y, XNew = NULL,
##                         nQ = 17L, nW = 32L, nR,
##                         deriv = 0L,
##                         checkX = TRUE) {

##     L3 <- min(c(n - 1, (3 * d + 1) / 2))
##     L5 <- min(c(n - 1, (5 * d + 1) / 2))
##     LM5 <- 3 * d + 2 * L

    

## }

## RIPPLE <- function(X, y, XNew = NULL,
##                    nQ = 17L, nW = 32L, nR,
##                    deriv = 0L,
##                    checkX = TRUE) {

##     L5 <- min(c(n - 1, (5 * d + 1) / 2))
##     LM5 <- 3 * d + 2 * L5

## }
