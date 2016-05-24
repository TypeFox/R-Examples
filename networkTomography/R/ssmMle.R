
#' Evaluate marginal log-likelihood for calibration SSM
#'
#' Evaluates marginal log-likelihood for calibration SSM of Blocker & Airoldi
#' (2011) using Kalman filtering. This is very fast and numerically stable,
#' using the univariate Kalman filtering and smoothing functions of \code{KFAS}
#' with Fortran implementations.
#'
#' @param theta numeric vector (length k+1) of parameters. theta[-1] =
#'      log(lambda), and theta[1] = log(phi)
#' @param Ft evolution matrix (k x k) for OD flows; include fixed
#       autoregressive parameters
#' @param yt matrix (k x n) of observed link loads, one observation per column
#' @param Zt observation matrix for system; should be routing matrix A
#' @param Rt covariance matrix for observation equation; typically small and
#'      fixed
#' @param k integer number of OD flows to infer
#' @param tau numeric power parameter for mean-variance relationship
#' @param initScale numeric inflation factor for time-zero state covariance;
#'      defaults to steady-state variance setting
#' @param nugget small positive value to add to diagonal of state evolution
#'      covariance matrix to ensure numerical stability
#' @return numeric marginal log-likelihood obtained via Kalman smoothing
#' @keywords models multivariate ts
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
#' time series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @export
#' @family calibrationModel
llCalibration <- function(theta, Ft, yt, Zt, Rt, k=ncol(Ft), tau=2,
                          initScale=1 / (1 - diag(Ft)^2),
                          nugget=sqrt(.Machine$double.eps)) {
    # Parse parameters
    lambda <- exp(theta[-1])
    phi <- exp(theta[1])

    P1 <- diag_mat(initScale*(phi*lambda^tau + nugget))
    a1 <- lambda / (1 - diag(Ft))

    # Setup matrices
    yt <- yt - as.vector(Zt %*% a1)
    yt <- as.ts(t(yt))
    V <- diag_mat(phi*lambda^tau + nugget)

    ssm <- SSModel(yt~SSMcustom(Z=Zt, T=Ft, R=diag_mat(rep(1, ncol(Ft))),
                                Q=V, a1=a1, P1=P1),
                   H=Rt, distribution='gaussian')
    return(logLik(ssm)) # S3 method for log-likelihood of SSModel object
}

#' Filtering & smoothing at MLE for calibration SSM
#'
#' Run Kalman filtering and smoothing at calculated MLE for parameters of
#' calibration SSM. This is used to obtain point and covariance estimates for
#' the actual OD flows X following estimation of other parameters.
#'
#' @param mle numeric vector (length k+1) of parameters. theta[-1] =
#'      log(lambda), and theta[1] = log(phi)
#' @param Ft evolution matrix (k x k) for OD flows; include fixed
#       autoregressive parameters
#' @param yt matrix (k x n) of observed link loads, one observation per column
#' @param Zt observation matrix for system; should be routing matrix A
#' @param Rt covariance matrix for observation equation; typically small and
#'      fixed
#' @param k integer number of OD flows to infer
#' @param tau numeric power parameter for mean-variance relationship
#' @param initScale numeric inflation factor for time-zero state covariance;
#'      defaults to steady-state variance setting
#' @param nugget small positive value to add to diagonal of state evolution
#'      covariance matrix to ensure numerical stability
#' @return numeric marginal log-likelihood obtained via Kalman smoothing
#' @return list containing result of Kalman smoothing; see \code{\link{SSModel}}
#'      and \code{\link{KFS}} for details
#' @keywords models multivariate ts
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
#' time series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @export
#' @family calibrationModel
mle_filter <- function(mle, Ft, yt, Zt, Rt,
                       k=ncol(Ft), tau=2, initScale=1 / (1 - diag(Ft)^2),
                       nugget=sqrt(.Machine$double.eps)) {
    # Parse parameters
    lambda <- exp(mle$par[-1])
    phi <- exp(mle$par[1])

    P1 <- diag_mat(initScale*(phi*lambda^tau + nugget))
    a1 <- lambda / (1 - diag(Ft))

    # Setup matrices
    yt <- yt - as.vector(Zt %*% a1)
    yt <- as.ts(t(yt))
    V <- diag_mat(phi*lambda^tau + nugget)

    # Run Kalman filter
    ssm <- SSModel(yt~SSMcustom(Z=Zt, T=Ft, R=diag_mat(rep(1, ncol(Ft))),
                                Q=V, a1=a1, P1=P1),
                   H=Rt, distribution='gaussian')
    filter.out <- KFS(ssm, smoothing='state')

    return(filter.out)
}

#' Estimation for the linear SSM calibration model of Blocker & Airoldi (2011)
#'
#' Maximum likelihood estimation of the parameters of the calibration model from
#' Blocker & Airoldi (2011) via direct numerical maximization of the marginal
#' log-likelihood. This relies upon efficient Kalman smoothing to evaluate the
#' marginal likelihood, which is provided here by the \code{KFAS} package.
#'
#' @param tme integer time at which to center moving window for estimation
#' @param y matrix (n x m) of observed link loads from all times (not just the
#'      window used for estimation; one observation per row
#' @param A routing matrix (m x k) for network; should be full row rank
#' @param Ft matrix (k x k) containing fixed autoregressive parameters for state
#'      evolution equation; upper-left block of overall matrix for expanded
#'      state
#' @param Rt covariance matrix for observation equation; typically small and
#'      fixed
#' @param lambda0 matrix (n x k) of initial estimates for lambda (e.g. obtained
#'      via IPFP)
#' @param phihat0 numeric vector (length n) of initial estimates for phi
#' @param tau numeric power parameter for mean-variance relationship
#' @param w number of observations to use for rolling-window estimation; handles
#'      boundary cases cleanly
#' @param initScale numeric inflation factor for time-zero state covariance;
#'      defaults to steady-state variance setting
#' @param nugget small positive value to add to diagonal of state evolution
#'      covariance matrix to ensure numerical stability
#' @param verbose logical to select verbose output from algorithm
#' @param logTrans logical whether to log-transform parameters for optimization.
#'      If FALSE, sets method to "L-BFGS-B".
#' @param method optimization method to use (in optim calls)
#' @param optimArgs list of arguments to append to control argument for optim.
#'      Can include all arguments except for fnscale, which is automatically set
#' @return list containing \code{lambdahat}, a numeric vector (length k)
#'      containing the MLE for lambda; \code{phihat}, the MLE for phi;
#'      \code{xhat}, the smoothed estimates of the OD flows for the window used
#'      as a k x w matrix; and \code{varhat}, a k x w matrix containing the 
#'      diagonal of the estimated covariance for each OD flow in the window
#' @keywords models multivariate ts
#' @references A.W. Blocker and E.M. Airoldi. Deconvolution of mixing
#' time series on a graph. Proceedings of the Twenty-Seventh Conference Annual
#' Conference on Uncertainty in Artificial Intelligence (UAI-11) 51-60, 2011.
#' @export
#' @family calibrationModel
#' @examples
#' data(bell.labs)
#' 
#' lambda0 <- matrix(1, nrow(bell.labs$Y), ncol(bell.labs$A))
#' lambda0[100,] <- ipfp(y=bell.labs$Y[100,], A=bell.labs$A,
#'                       x0=rep(1, ncol(bell.labs$A)))
#' phihat0 <- rep(1, nrow(bell.labs$Y))
#' Ft <- 0.5 * diag_mat(rep(1, ncol(bell.labs$A)))
#' Rt <- 0.01 * diag_mat(rep(1, nrow(bell.labs$A)))
#' 
#' # Not run
#' #fit.calibration <- calibration_ssm(tme=100, y=bell.labs$Y, A=bell.labs$A,
#' #                                   Ft=Ft, Rt=Rt, lambda0=lambda0,
#' #                                   phihat0=phihat0, w=23)
calibration_ssm <- function(tme, y, A, Ft, Rt, lambda0, phihat0, tau=2, w=11,
                            initScale=1 / (1 - diag(Ft)^2),
                            nugget=sqrt(.Machine$double.eps), verbose=FALSE,
                            logTrans=TRUE, method="L-BFGS-B",
                            optimArgs=list()) {
    # Calculate dimensions
    k <- ncol(A)
    l <- ncol(y)

    # Calculate window parameters
    h <- floor(w/2)

    # Calculate length of window and index of tme within window
    if ( (tme-h>0) & (tme+h<nrow(y)) ) {
        n <- w
        t_ind <- h+1
    } else {
        # Handle border cases
        if (tme-h<=0) {
            n <- tme+h
            t_ind <- tme
        } else {
            n <- nrow(y)-tme+h+1
            t_ind <- h+1
        }
    }

    # Print tme and window length if verbose
    if (verbose) {
        cat(sprintf('tme = %d; n = %d\n', tme, n), file=stderr())
    }

    # Setup data
    yt <- t( y[max(1,tme-h):min(tme+h,nrow(y)),] )
    lambda <- lambda0[tme,]
    phi <- phihat0[tme]

    # Print starting values if verbose
    if (verbose) {
        cat('Starting value for lambda:\n', file=stderr())
        print(lambda)
        cat(sprintf('Starting value for phi:\t%g\n', phi), file=stderr())
    }

    # Run numerical optimization using prediction error formulation
    # of Airoldi 2003 (section 3.1.2); with a Fortran implementation
    # of the univariate Kalman filter and smoother of Koopman & Durbin (2000,
    # 2003), this is far more efficient than the EM iterations (quadratic vs.
    # linear convergence, plus much less memory usage in the smoothing stage).
    if (logTrans) {
        obj <- llCalibration
        theta0 <- log(c(phi, lambda))
        lower <- -Inf
        upper <- Inf
    } else {
        obj <- function(theta, ...) llCalibration(log(theta), ...)
        theta0 <- c(phi, lambda)
        lower <- c(sqrt(.Machine$double.eps), rep(1, k))
        upper <- max(yt)
        method <- "L-BFGS-B"
    }

    # Print initial value for objective function if verbose
    if (verbose)
        cat(sprintf('Initial value of objective function:\t%g\n', 
                    obj(theta0, Ft=Ft, yt=yt, Zt=A, Rt=Rt, tau=tau,
                        initScale=initScale, nugget=nugget)), file=stderr())

    mle <- optim(par=theta0, fn=obj,
                 Ft=Ft, yt=yt, Zt=A, Rt=Rt, tau=tau,
                 initScale=initScale,
                 nugget=nugget,
                 method=method,
                 lower=lower, upper=upper,
                 control=c(list(fnscale=-1), optimArgs))
    
    if (!logTrans)
        mle$par <- log(mle$par)

    # Print optim diagnostics if verbose
    if (verbose) {
        cat(sprintf("Method: %s\n", method), file=stderr())
        cat(sprintf('Convergence code: %d\n', mle$convergence), file=stderr())
        cat('Function evaluations:\n', file=stderr())
        print(mle$counts)
    }

    # Obtain Kalman filter output at MLE
    f.out <- mle_filter(mle=mle, Ft=Ft, yt=yt, Zt=A, Rt=Rt, tau=tau,
                        initScale=initScale, nugget=nugget)
    varhat <- apply(f.out$P, 3, diag)

    # Reformat for output
    lambdahat   <- exp(mle$par[-1])
    xhat        <- f.out$alphahat[1:k,t_ind] + lambdahat / (1 - diag(Ft))
    phihat      <- exp(mle$par[1])
    varhat      <- varhat[1:k,t_ind]

    # Return results
    return(list(lambdahat=lambdahat, phihat=phihat, xhat=xhat, varhat=varhat))
}

