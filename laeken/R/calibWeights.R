# ---------------------------------------
# Author: Andreas Alfons
#         Vienna University of Technology
# ---------------------------------------

#' Calibrate sample weights
#' 
#' Calibrate sample weights according to known marginal population totals.
#' Based on initial sample weights, the so-called \emph{g}-weights are computed
#' by generalized raking procedures.
#' 
#' The final sample weights need to be computed by multiplying the resulting
#' \emph{g}-weights with the initial sample weights.
#' 
#' @encoding utf8
#' 
#' @param X a matrix of binary calibration variables (see
#' \code{\link{calibVars}}).
#' @param d a numeric vector giving the initial sample weights.
#' @param totals a numeric vector of population totals corresponding to the
#' calibration variables in \code{X}.
#' @param q a numeric vector of positive values accounting for
#' heteroscedasticity.  Small values reduce the variation of the
#' \emph{g}-weights.
#' @param method a character string specifying the calibration method to be
#' used.  Possible values are \code{"linear"} for the linear method,
#' \code{"raking"} for the multiplicative method known as raking and
#' \code{"logit"} for the logit method.
#' @param bounds a numeric vector of length two giving bounds for the g-weights
#' to be used in the logit method.  The first value gives the lower bound (which
#' must be smaller than or equal to 1) and the second value gives the upper
#' bound (which must be larger than or equal to 1).
#' @param maxit a numeric value giving the maximum number of iterations.
#' @param tol the desired accuracy for the iterative procedure.
#' @param eps the desired accuracy for computing the Moore-Penrose generalized
#' inverse (see \code{\link[MASS]{ginv}}).
#' 
#' @return A numeric vector containing the \emph{g}-weights.
#' 
#' @note This is a faster implementation of parts of
#' \code{\link[sampling]{calib}} from package \code{sampling}.  Note that the
#' default calibration method is raking and that the truncated linear method is
#' not yet implemented.
#' 
#' @author Andreas Alfons
#' 
#' @seealso \code{\link{calibVars}}, \code{\link{bootVar}}
#' 
#' @references Deville, J.-C. and \enc{Särndal}{Saerndal}, C.-E. (1992)
#' Calibration estimators in survey sampling. \emph{Journal of the American
#' Statistical Association}, \bold{87}(418), 376--382.
#' 
#' Deville, J.-C., \enc{Särndal}{Saerndal}, C.-E. and Sautory, O. (1993)
#' Generalized raking procedures in survey sampling. \emph{Journal of the
#' American Statistical Association}, \bold{88}(423), 1013--1020.
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' # construct auxiliary 0/1 variables for genders
#' aux <- calibVars(eusilc$rb090)
#' # population totals
#' totals <- c(3990798, 4191431)
#' # compute g-weights
#' g <- calibWeights(aux, eusilc$rb050, totals)
#' # compute final weights
#' weights <- g * eusilc$rb050
#' summary(weights)
#' 
#' @export
#' @import MASS

calibWeights <- function(X, d, totals, q = NULL, 
        method = c("raking", "linear", "logit"), 
        bounds = c(0, 10), maxit = 500, tol = 1e-06, 
        eps = .Machine$double.eps) {
    
    ## initializations and error handling
    X <- as.matrix(X)
    d <- as.numeric(d)
    totals <- as.numeric(totals)
    haveNA <- c(any(is.na(X)), any(is.na(d)), 
        any(is.na(totals)), !is.null(q) && any(is.na(q)))
    if(any(haveNA)) {
        argsNA <- c("'X'", "'d'", "'totals'", "'q'")[haveNA]
        stop("missing values in the following arguments", 
            paste(argsNA, collapse=", "))
    }
    n <- nrow(X)  # number of rows
    if(length(d) != n) stop("length of 'd' not equal to number of rows in 'X'")
    p <- ncol(X)  # number of columns
    if(length(totals) != p) {
        stop("length of 'totals' not equal to number of columns in 'X'")
    }
    if(is.null(q)) q <- rep.int(1, n)
    else {
        q <- as.numeric(q)
        if(length(q) != n) {
            stop("length of 'q' not equal to number of rows in 'X'")
        }
        if(any(is.infinite(q))) stop("infinite values in 'q'")
    }
    method <- match.arg(method)
    
    ## computation of g-weights
    if(method == "linear") {
        ## linear method (no iteration!)
        lambda <- ginv(t(X * d * q) %*% X, tol=eps) %*% (totals - as.vector(t(d) %*% X))
        g <- 1 + q * as.vector(X %*% lambda)  # g-weights
    } else {
        ## multiplicative method (raking) or logit method
        lambda <- matrix(0, nrow=p)  # initial values
        # function to determine whether teh desired accuracy has 
        # not yet been reached (to be used in the 'while' loop)
        tolNotReached <- function(X, w, totals, tol) {
            max(abs(crossprod(X, w) - totals)/totals) >= tol
        }
        if(method == "raking") {
            ## multiplicative method (raking)
            # some initial values
            g <- rep.int(1, n)  # g-weights
            w <- d  # sample weights
            ## iterations
            i <- 1
            while(!any(is.na(g)) && tolNotReached(X, w, totals, tol) && i <= maxit) {
                # here 'phi' describes more than the phi function in Deville, 
                # Saerndal and Sautory (1993); it is the whole last term of 
                # equation (11.1)
                phi <- t(X) %*% w - totals
                T <- t(X * w)
                dphi <- T %*% X  # derivative of phi function (to be inverted)
                lambda <- lambda - ginv(dphi, tol=eps) %*% phi  # update 'lambda'
                g <- exp(as.vector(X %*% lambda) * q)  # update g-weights
                w <- g * d  # update sample weights
                i <- i + 1  # increase iterator
            }
            ## check wether procedure converged
            if(any(is.na(g)) || i > maxit) {
                warning("no convergence")
                g <- NULL
            }
        } else {
            ## logit (L, U) method
            ## error handling for bounds
            if(length(bounds) < 2) stop("'bounds' must be a vector of length 2")
            else bounds <- bounds[1:2]
            if(bounds[1] >= 1) stop("the lower bound must be smaller than 1")
            if(bounds[2] <= 1) stop("the lower bound must be larger than 1")
            ## some preparations
            A <- diff(bounds)/((1 - bounds[1]) * (bounds[2] - 1))
            # function to bound g-weights
            getG <- function(u, bounds) {
                (bounds[1] * (bounds[2]-1) + bounds[2] * (1-bounds[1]) * u) / 
                    (bounds[2]-1 + (1-bounds[1]) * u)
            }
            ## some initial values
            g <- getG(rep.int(1, n), bounds)  # g-weights
            # in the procedure, g-weights outside the bounds are moved to the 
            # bounds and only the g-weights within the bounds are adjusted.
            # these duplicates are needed since in general they are changed in 
            # each iteration while the original values are also needed
            X1 <- X
            d1 <- d
            totals1 <- totals
            q1 <- q
            g1 <- g
            indices <- 1:n
            # function to determine which g-weights are outside the bounds
            anyOutOfBounds <- function(g, bounds) {
                any(g < bounds[1]) || any(g > bounds[2])
            }
            ## iterations
            i <- 1
            while(!any(is.na(g)) && (tolNotReached(X, g*d, totals, tol) ||
                    anyOutOfBounds(g, bounds)) && i <= maxit) {
                # if some of the g-weights are outside the bounds, these values 
                # are moved to the bounds and only the g-weights within the 
                # bounds are adjusted
                if(anyOutOfBounds(g, bounds)) {
                    g[g < bounds[1]] <- bounds[1]
                    g[g > bounds[2]] <- bounds[2]
                    # values within the bounds
                    tmp <- which(g > bounds[1] & g < bounds[2])
                    if(length(tmp) > 0) {
                        indices <- tmp
                        X1 <- X[indices,]
                        d1 <- d[indices]
                        if(length(indices) < n) {
                            totals1 <- totals - as.vector(t(g[-indices] * d[-indices]) %*% X[-indices, , drop=FALSE])
                        }
                        q1 <- q[indices]
                        g1 <- g[indices]
                    }
                }
                w1 <- g1 * d1  # current sample weights
                # here 'phi' describes more than the phi function in Deville, 
                # Saerndal and Sautory (1993); it is the whole last term of 
                # equation (11.1)
                phi <- t(X1) %*% w1 - totals1
                T <- t(X1 * w1)
                dphi <- T %*% X1  # derivative of phi function (to be inverted)
                lambda <- lambda - ginv(dphi, tol=eps) %*% phi  # update 'lambda'
                # update g-weights
                u <- exp(A * as.vector(X1 %*% lambda) * q1)
                g1 <- getG(u, bounds)
                g[indices] <- g1
                i <- i+1  # increase iterator
            }
            ## check wether procedure converged
            if(any(is.na(g)) || i > maxit) {
                warning("no convergence")
                g <- NULL
            }
        }
        
    } 
    
    ## return g-weights
    return(g)
}
