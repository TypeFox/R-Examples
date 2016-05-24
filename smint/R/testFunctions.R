##' Branin-Hoo 2-dimensional test function.
##'
##' The Branin-Hoo function is defined here over \eqn{[0,\,1] \times
##' [0,\,1]}{[0, 1] x [0, 1]}, instead of \eqn{[-5,\,0] \times
##' [10,\,15]}{[-5, 0] x [10, 15]} as usual.  It has 3 global minima
##' at (nearly) : \eqn{\mathbf{x}_1 = [0.96, \,0.15]^\top}{x1 = c(0.96, 0.15)},
##' \eqn{\mathbf{x}_2 = [0.12, \,0.82]^\top}{x2 = c(0.12, 0.82)}
##' and \eqn{\mathbf{x}_3 = [0.54,\,0.15]^\top}{x3 = c(0.54, 0.15)}.
##'
##' @title  Branin-Hoo 2-dimensional test function
##'
##' @param x Numeric vector with length 2.
##' 
##' @return The Branin-Hoo function's value.
##'
##' @author David Ginsbourger
##'
##' @examples
##' GD <- Grid(nlevels = c("x" = 20, "y" = 20))
##' x <- levels(GD)[[1]]; y <- levels(GD)[[2]]
##' f  <- apply_Grid(GD, branin)
##' dim(f) <- nlevels(GD)
##' contour(x = x, y = y, z = f, nlevels = 40)
##' nOut <- 100; Xout2 <- array(runif(nOut * 2), dim = c(nOut, 2))
##' colnames(Xout2) <- c("x", "y")
##' 
##' ## interpolate using default method (Lagrange) 
##' GIL <- interp_Grid(X = GD, Y = f, Xout = Xout2)
##'
##' ## interpolate using a natural spline
##' GIS <- interp_Grid(X = GD, Y = f, Xout = Xout2,
##'                  cardinalBasis1d = function(x, xout) {
##'                      cardinalBasis_natSpline(x = x, xout = xout)$CB
##'                   })
##' F <- apply(Xout2, 1, branin)
##' mat <- cbind(Xout2, fTrue = F, fIntL = GIL, errorLag = F - GIL,
##'              fIntS = GIS, errorSpline = F - GIS)
##' apply(mat[ , c("errorLag", "errorSpline")], 2, function(x) mean(abs(x)))
##' \dontrun{
##' ## for the users of the "rgl" package only...
##' library(rgl)
##' persp3d(x = x, y = y, z = f, aspect = c(1, 1, 0.5), col = "lightblue", alpha = 0.8)
##' spheres3d(Xout2[ , 1], Xout2[ , 2], GIS, col = "orangered",
##'           radius = 2)
##' }
branin <- function (x) {
    x1 <- x[1] * 15 - 5
    x2 <- x[2] * 15
    (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 
        1/(8 * pi)) * cos(x1) + 10
}

##' Test functions for/from SHEPPACK
##' 
##' These functions are described in the article cited in the \bold{references}
##' section.
##' \deqn{f_1(\mathbf{x}) = 1 + \frac{2}{d} \, \left| d/2 - \sum_i x_i
##' \right|}{ f1(x) = 1 + (2/d) | (d/2) - sum_i x_i |}
##' \deqn{f_2(\mathbf{x}) = 1 - \frac{2}{d}\, \sum_i \left|x_i - 0.5
##' \right|}{ f2(x) = 1 - (2/d) sum_i |x_i - 0.5 |}
##' \deqn{f_3(\mathbf{x}) = 1 - 2 \, \max_{i} \left|x_i - 0.5 \right|}{
##' f3(x) = 1 - 2 max_{i} |x_i - 0.5 |}
##' \deqn{f_4(\mathbf{x}) = \prod_i \left[ 1 - 2 \left| x_i - 0.5
##' \right| \right]}{ f4(x) = prod_i [ 1 - 2 | x_i - 0.5 | ]}
##' \deqn{f_5(\mathbf{x}) = 1 - c_5 \, \left[ \sum_i \left|x_i - 0.5
##' \right| + \prod_i \left| x_i - 0.5 \right| \right]}{ f5(x) = 1 -
##' c_5 [ sum_i |x_i - 0.5 | + prod_i | x_i - 0.5 | ]}
##' where \eqn{c_5 = d/2 + (0.5)^d}, and all sums or products are for
##' \eqn{i=1} to \eqn{d}. All these functions are defined on
##' \eqn{[0,\,1]^d} and take values in \eqn{[0,1]}. The four functions
##' \eqn{f_i} for \eqn{i > 1} have an unique maximum located at
##' \eqn{\mathbf{x}^\star}{xStar} with all coordinates \eqn{x_j^\star
##' = 0.5}{xStar[j] = 0.5} and \eqn{f_i(\mathbf{x}^\star) =1}{f_i(xStar) = 1}.
##' 
##' @aliases ShepFun2 ShepFun3 ShepFun4 ShepFun5 ShepFuns
##' 
##' @title Test functions for/from SHEPPACK
##'
##' @param x A numeric vector with arbitrary length. 
##'
##' @return Function's value.
##'
##' @references W.I. Thacker, J. Zhang, L.T. Watson, J.B. Birch,
##' M.A. Iyer and M.W. Berry (2010). Algorithm 905: SHEPPACK: Modified
##' Shepard Algorithm for Interpolation of Scattered Multivariate Data
##' \emph{ACM Trans. on Math. Software} (TOMS) Vol. 37, n. 3.
##' \href{http://dl.acm.org/citation.cfm?id=1824812}{link}
##' 
##' @note These functions are also exported as elements of the
##' \code{ShepFuns} list.
##' 
##' @examples
##' ## interpolate 'Shepfun3' for d = 4
##' d <- 4
##' GDd <- Grid(nlevels = rep(8, d))
##' fGrid <- apply_Grid(GDd, ShepFun3)
##' Xoutd <- matrix(runif(200 * d), ncol = d)
##' GI <- interp_Grid(X = GDd, Y = fGrid, Xout = Xoutd)
##' F <- apply(Xoutd, 1, ShepFun3)
##' max(abs(F - GI))
##'
##' ## 3D plot
##' require(lattice) 
##' X <- as.data.frame(Grid(nlevels = c("x1" = 30, "x2" = 30)))
##' df <- data.frame(x1 = numeric(0), x2 = numeric(0),
##'                  f = numeric(0), i = numeric(0))
##' for (i in 1:5) {
##'    f <- apply(X, 1, ShepFuns[[i]])
##'    df <- rbind(df, data.frame(x1 = X$x1, x2 = X$x2, f = f, i = i))
##' }
##' pl <- wireframe(f ~ x1 * x2 | i, data = df,
##'                 outer = TRUE, shade = FALSE, zlab = "",
##'                 screen = list(z = 20, x = -30),
##'                 strip = strip.custom(strip.names = c(TRUE),
##'                                      strip.levels = c(TRUE)),
##'                 main = "", horizontal = TRUE, col = "SpringGreen4")
##' pl

ShepFun1 <- function(x) {
    d <- length(x)
    f <- 2 * sum(x) / d
        if (f <= 1) return(f)
        else return(2 - f)
}

ShepFun2 <- function(x) {
    d <- length(x)
    f <- 1 - 2 * sum(abs(x - 0.5)) / d
}

ShepFun3 <- function(x) {
    f <- 1 - 2 * max(abs(x - 0.5)) 
}

ShepFun4 <- function(x) {
    prod(ifelse(x <= 0.5, 2*x, 2 * (1-x)))
}

ShepFun5 <- function(x) {
    d <- length(x) 
    xf <- abs(x - 0.5)
    1 - ( sum(xf) + prod(xf) ) / (d * 0.5 + 0.5^d)
}

ShepFuns <- list(ShepFun1, ShepFun2, ShepFun3, ShepFun4, ShepFun5)
