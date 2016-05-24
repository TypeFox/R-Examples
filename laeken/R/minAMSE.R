# ----------------------------------------
# Authors: Josef Holzer and Andreas Alfons
#          Vienna University of Technology
# ----------------------------------------

## nonlinear integer minimization is done by brute force
## it is strongly recommended to set bounds 'kmax' and 'mmax'

#' Weighted asymptotic mean squared error (AMSE) estimator
#' 
#' Estimate the scale and shape parameters of a Pareto distribution with an
#' iterative procedure based on minimizing the weighted asymptotic mean squared
#' error (AMSE) of the Hill estimator.
#' 
#' The weights used in the weighted AMSE depend on a nuisance parameter
#' \eqn{\rho}{rho}.  Both the optimal number of observations in the tail and the
#' nuisance parameter \eqn{\rho}{rho} are estimated iteratively using nonlinear
#' integer minimization.  This is currently done by a brute force algorithm,
#' hence it is stronly recommended to supply upper bounds \code{kmax} and
#' \code{mmax}.
#' 
#' See the references for more details on the iterative algorithm.
#' 
#' @aliases print.minAMSE
#' 
#' @param x for \code{minAMSE}, a numeric vector.  The \code{print} method is
#' called by the generic function if an object of class \code{"minAMSE"} is
#' supplied.
#' @param weight a character vector specifying the weighting scheme to be used
#' in the procedure.  If \code{"Bernoulli"}, the weight functions as described
#' in the \emph{Bernoulli} paper are applied.  If \code{"JASA"}, the weight
#' functions as described in the \emph{Journal of the Americal Statistical
#' Association} are used.
#' @param kmin An optional integer giving the lower bound for finding the
#' optimal number of observations in the tail.  It defaults to
#' \eqn{[\frac{n}{100}]}{[n/100]}, where \eqn{n} denotes the number of
#' observations in \code{x} (see the references).
#' @param kmax An optional integer giving the upper bound for finding the
#' optimal number of observations in the tail (see \dQuote{Details}).
#' @param mmax An optional integer giving the upper bound for finding the
#' optimal number of observations for computing the nuisance parameter
#' \eqn{\rho}{rho} (see \dQuote{Details} and the references).
#' @param tol an integer giving the desired tolerance level for finding the
#' optimal number of observations in the tail.
#' @param maxit a positive integer giving the maximum number of iterations.
#' @param \dots additional arguments to be passed to
#' \code{\link[base]{print.default}}.
#' 
#' @return An object of class \code{"minAMSE"} containing the following
#' components:
#' @returnItem kopt the optimal number of observations in the tail.
#' @returnItem x0 the corresponding threshold.
#' @returnItem theta the estimated shape parameter of the Pareto distribution.
#' @returnItem MSEmin the minimal MSE.
#' @returnItem rho the estimated nuisance parameter.
#' @returnItem k the examined range for the number of observations in the tail.
#' @returnItem MSE the corresponding MSEs.
#' 
#' @author Josef Holzer and Andreas Alfons
#' 
#' @seealso \code{\link{thetaHill}}
#' 
#' @references Beirlant, J., Vynckier, P. and Teugels, J.L. (1996) Tail index
#' estimation, Pareto quantile plots, and regression diagnostics. \emph{Journal
#' of the American Statistical Association}, \bold{91}(436), 1659--1667.
#' 
#' Beirlant, J., Vynckier, P. and Teugels, J.L. (1996) Excess functions and
#' estimation of the extreme-value index. \emph{Bernoulli}, \bold{2}(4),
#' 293--318.
#' 
#' Dupuis, D.J. and Victoria-Feser, M.-P. (2006) A robust prediction error
#' criterion for Pareto modelling of upper tails. \emph{The Canadian Journal of
#' Statistics}, \bold{34}(4), 639--658.
#' 
#' @keywords manip
#' 
#' @examples
#' data(eusilc)
#' # equivalized disposable income is equal for each household
#' # member, therefore only one household member is taken
#' minAMSE(eusilc$eqIncome[!duplicated(eusilc$db030)], 
#'     kmin = 50, kmax = 150, mmax = 250)
#' 
#' @export

minAMSE <- function(x, weight = c("Bernoulli", "JASA"), 
        kmin, kmax, mmax, tol = 0, maxit = 100) {
    ## initializations
    if(!is.numeric(x) || length(x) == 0) stop("'x' must be a numeric vector")
    if(any(i <- is.na(x))) x <- x[!i]
    x <- sort(x)
    n <- length(x)
    if(n == 0) stop("no observed values")
    weight <- match.arg(weight)
    kbounds <- c(trunc(n/100), n-2)
    mbound <- n-1
    if(missing(kmin)) kmin <- kbounds[1]
    if(missing(kmax)) kmax <- kbounds[2]
    if(missing(mmax)) mmax <- mbound
    if(!is.numeric(kmin) || length(kmin) == 0 || kmin[1] < 1) {
        stop("'kmin' must be a single positive integer")
    } else kmin <- kmin[1]
    if(!is.numeric(kmax) || length(kmax) == 0 || kmax[1] <= kmin) {
        stop("'kmax' must be a single positive integer larger than 'kmin'")
    } else kmax <- kmax[1]
    if(!is.numeric(mmax) || length(mmax) == 0 || mmax[1] <= kmax) {
        stop("'mmax' must be a single positive integer larger than 'kmax'")
    } else mmax <- mmax[1]
    if(!is.numeric(maxit) || length(maxit) == 0 || maxit[1] < 1) {
        stop("'maxit' must be a single positive integer")
    } else maxit <- maxit[1]
    ## check bounds for k
    if(kmin < kbounds[1]) {
        kmin <- kbounds[1]
        warning("'kmin' is set to ", kbounds[1], 
            ", as this is the suggested minumum")
    }
    if(kmax > kbounds[2]) {
        kmax <- kbounds[2]
        warning("'kmax' is set to ", kbounds[2], 
            ", as this is the allowed maximum")
    }
    ## check bound for m
    if(mmax > mbound) {
        mmax <- mbound
        warning("'mmax' is set to ", mbound, 
            ", as this is the allowed maximum")
    }
    ## Hill estimates of theta for range of k 
    kl <- trunc(kmin/2)
    theta <- rep.int(NA, kmax)
    theta[kl:mmax] <- sapply(kl:mmax, function(k) thetaHill(x, k))  # shape
    ## initial estimate of k
    k <- kmin:kmax  # range of k to search for minimum
    MSE <- mapply(function(k, theta) MSEinit(x, k, theta), 
        k, theta[(kmin:kmax)-kl+1])
    k0 <- k[which.min(MSE)]
    theta0 <- theta[k0]
    ## initial estimate of rho
    m <- (k0+1):mmax
    rho <- sapply(m, function(m) Rm(x, theta, m, k0))
    cr <- mapply(function(m, rho) critRm(x, rho, theta, m, k0), m, rho)
    rho0 <- rho[which.min(cr)]
    ## iterative procedure
    for(i in 1:maxit) {
        ## estimate k
        MSE <- sapply(k, function(k) MSEopt(x, k, theta0, rho0, weight))
        tmp <- which.min(MSE)
        MSEmin <- MSE[tmp]
        kn <- k[tmp]
        thetan <- theta[kn]
        ## estimate rho
        m <- (kn+1):mmax
        rho <- sapply(m, function(m) Rm(x, theta, m, kn))
        cr <- mapply(function(m, rho) critRm(x, rho, theta, m, k0), m, rho)
        rhon <- rho[which.min(cr)]
        if(abs(kn-k0) <= tol) break
        else {
            k0 <- kn
            theta0 <- thetan
            rho0 <- rhon
        }
    }
    ## return results
    res <- list(kopt=kn, x0=x[n-kn], theta=thetan, 
        MSEmin=MSEmin, rho=rhon, k=k, MSE=MSE)
    class(res) <- "minAMSE"
    res
}


## internal functions for the evaluation of the MSEopt criterion
## x is not expected to contain missing values and is assumed to be sorted

MSEinit <- function(x, k, theta) {
    n <- length(x)
    x0 <- x[n-k]  # threshold (scale parameter)
    y <- log(x[(n-k+1):n]/x0)  # relative excesses
    nyhat <- log((k:1)/(k+1))/theta  # negative predicted values
    ## MSE
    1/k * sum((y + nyhat)^2)
}

MSEopt <- function(x, k, theta, rho, weight = c("Bernoulli", "JASA")) {
    n <- length(x)
    x0 <- x[n-k]  # threshold (scale parameter)
    y <- log(x[(n-k+1):n]/x0)  # relative excesses
    h <- k:1
    nyhat <- log(h/(k+1))/theta  # negative predicted values
    ## weight functions according to paper in Bernoulli or JASA
    i <- 1:k
    hv <- i/(k+1)
    if(weight == "Bernoulli") {
        wk1 <- hv
        wk2 <- -log(i/(k+1))
    } else {
        wk1 <- rep.int(1, k)
        wk2 <- h/(k+1)  # second weight function (first is identical to 1)
    }
    ## define delta functions
    tmp1 <- hv^(-1)-1
    tmp2 <- (1-rho)^2
    tmp3 <- ((hv^(-rho)-1)/rho)^2
    ak1 <- mean(wk1*tmp1)
    ak2 <- mean(wk2*tmp1)
    bk1 <- tmp2 * mean(wk1*tmp3)
    bk2 <- tmp2 * mean(wk2*tmp3)
    den <- (ak1*bk2 - bk1*ak2)  # denomitator for delta functions
    delta1 <- (bk2 - ak2) / den
    delta2 <- (ak1 - bk1) / den
    ## define optimal weight function
    woptk <- delta1*wk1 + delta2*wk2
    ## WMSE
    mean(woptk * (y + nyhat)^2)
}


## internal functions for estimating rho
## requirements for m and k are assumed to be fulfilled

Rm <- function(x, theta, m, k) {
    mk <- m+k
    # denominators 2 and 4 do not cause problems with floating point arithmetic
    Hmk4 <- theta[trunc(mk/4)]
    Hmk2 <- theta[trunc(mk/2)]
    Hm2 <- theta[trunc(m/2)]
    Hm <- theta[m]
    log(abs((Hmk4-Hmk2)/(Hm2-Hm))) / (log(2*m/(m-k)))
}


## internal functions for the evaluation of the criterion for rho
## x is not expected to contain missing values and is assumed to be sorted

critRm <- function(x, rho, theta, m, k) {
    j <- k:(m-1)
    l <- log(abs((theta[trunc(j/2)] - theta[j])/(theta[trunc(m/2)] - theta[m])))
    mean((l - rho*log(m/j))^2)
}
