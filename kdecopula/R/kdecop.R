#' Bivariate kernel copula density estimation
#' 
#' Based on samples from a bivariate copula, the copula density is estimated.
#' The user can choose between different methods. If no bandwidth is provided
#' by the user, it will be set by a method-specific automatic selection
#' procedure. The related (d/p/r)kdecop functions evaluate the density and cdf 
#' or simulate synthetic data, respectively.
#'   
#' @param udata \code{nx2} matrix of copula data.
#' @param bw bandwidth specification; if \code{NA}, \code{bw} is selected
#' automatically; Otherwise, please provide (for the respective method) \cr
#' \code{"MR", "beta"}: a positive real number, \cr
#' \code{"T"}: a \eqn{2x2} matrix for method, \cr
#' \code{"TLL1", "TLL2"}: a list with (named) entries \code{B} and \code{alpha}
#' containing the \eqn{2x2} rotation matrix \code{B} and the nearest-neighbor 
#' fraction \code{alpha}, \cr
#' \code{"TTCV", "TTPI"}: a numeric vector of length four containing \eqn{(h, 
#' \rho, \theta_1, \theta_2)}, c.f. Wen and Wu (2015).
#' @param mult bandwidth multiplier, has to be positive; useful for making 
#' estimates more/less smooth manually.
#' @param method 
#' \code{"MR"}: mirror-reflection estimator, \cr 
#' \code{"beta"}: beta kernel estimator, \cr 
#' \code{"T"}: transformation estimator with standard bivariate kernel 
#' estimation, \cr 
#' \code{"TLL1"}: transformation estimator with log-linear local likelihood
#' estimation and nearest-neighbor bandwidths (Geenenens et al., 2014), \cr 
#' \code{"TLL2"}: transformation estimator with log-quadradtic local likelihood 
#' estimation and nearest-neighbor bandwidths (Geenenens et al., 2014), \cr
#' \code{"TTPI"}: tapered transformation estimator with plug-in bandwidths, \cr
#' \code{"TTCV"}: tapered transformation estimator with profile cross-validation
#' bandwidths.
#' @param knots integer; number of knots in each dimension for the spline
#' approximation.
#' @param renorm.iter integer; number of iterations for the renormalization
#' procedure (see \emph{Details}).
#' @param info logical; if \code{TRUE}, additional information about the
#' estimate will be gathered (see \emph{Value}).
#' 
#' 
#' @return The function \code{\link[kdecopula:kdecop]{kdecop}} returns an
#' object of class \code{kdecopula} that contains all information necessary for
#' evaluation of the estimator. If no bandwidth was provided in the function
#' call, the automatically selected value can be found in the variable
#' \code{object$bw}. If \code{info=TRUE}, also the following will be available
#' under \code{object$info}: 
#' \item{likvalues}{Estimator evaluated in sample points} 
#' \item{loglik}{Log likelihood} 
#' \item{effp}{Effective number of parameters} 
#' \item{AIC}{Akaike information criterion}
#' \item{cAIC}{Bias-corrected version of Akaike information criterion}
#' \item{BIC}{Bayesian information criterion.} \cr 
#' The density estimate can be evaluated on arbitrary points with 
#' \code{\link[kdecopula:dkdecop]{dkdecop}}; the cdf with 
#' \code{\link[kdecopula:pkdecop]{pkdecop}}. Furthermore, synthetic data can be
#' simulated with \code{\link[kdecopula:rkdecop]{rkdecop}}, and several plotting
#' options are available with \code{\link[kdecopula:plot.kdecopula]{plot}}
#' and \code{\link[kdecopula:contour.kdecopula]{contour}}.
#' 
#' @details Details on the estimation methods and bandwidth selection can be
#' found in Geenens et al. (2014) for methods \code{TLL1/2} and Nagler (2014) 
#' for other methods. We use a Gaussian product kernel function for all methods 
#' except the beta kernel estimator.\cr 
#' 
#' Kernel estimates are usually no proper copula densities. In particular, the
#' estimated marginal densities are not uniform. We mitigate this issue bei
#' implementing a renormalization procedure. The number of iterations of the
#' renormalization algorithm can be specified with the \code{renorm.iter}
#' argument. Typically, a very small number of iterations is sufficient. \cr
#' 
#' The implementation of the tapered transformation estimator ("TTPI"/"TTCV") 
#' was kindly provided by Kuangyu Wen. 
#' 
#' @author Thomas Nagler
#' 
#' @seealso 
#' \code{\link[kdecopula:kdecopula]{kdecopula}},
#' \code{\link[kdecopula:plot.kdecopula]{plot.kdecopula}},
#' \code{\link[kdecopula:dkdecop]{dkdecop}},
#' \code{\link[kdecopula:pkdecop]{pkdecop}},
#' \code{\link[kdecopula:rkdecop]{rkdecop}}
#' 
#' @references 
#' Geenens, G., Charpentier, A., and Paindaveine, D. (2014).
#' Probit transformation for nonparametric kernel estimation of the copula
#' density.
#' arXiv:1404.4414 [stat.ME]. 
#' \cr \cr 
#' Nagler, T. (2014). 
#' Kernel Methods for Vine Copula Estimation.
#' Master's Thesis, Technische Universitaet Muenchen,
#' \url{https://mediatum.ub.tum.de/node?id=1231221} 
#' \cr \cr
#' Wen, K. and Wu, X. (2015).
#' Transformation-Kernel Estimation of the Copula Density,
#' Working paper,
#' \url{http://agecon2.tamu.edu/people/faculty/wu-ximing/agecon2/public/copula.pdf}
#' 
#' @examples
#' 
#' ## load data and transform with empirical cdf
#' data(wdbc)
#' udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))
#' 
#' ## estimation of copula density of variables 5 and 6
#' dens.est <- kdecop(udat[, 5:6])
#' summary(dens.est)
#' plot(dens.est) 
#' 
#' ## evaluate density estimate at (u1,u2)=(0.123,0.321)
#' dkdecop(c(0.123, 0.321), dens.est) 
#' 
#' ## evaluate cdf estimate at (u1,u2)=(0.123,0.321)
#' pkdecop(c(0.123, 0.321), dens.est) 
#' 
#' ## simulate 500 samples from density estimate
#' plot(rkdecop(500, dens.est))  # pseudo-random
#' plot(rkdecop(500, dens.est), quasi = TRUE)  # quasi-random
#' 
kdecop <- function(udata, bw = NA, mult = 1, method = "TLL2", knots = 50, renorm.iter = 3L, info = TRUE) {
    udata <- as.matrix(udata)
    n <- nrow(udata)
    d <- ncol(udata)
    
    ## sanity checks
    if (n < 2)
        stop("Number of observations has to be at least 2.")
    if (d != 2)
        stop("Dimension has to be 2.")
    if (any(udata > 1) | any(udata < 0))
        stop("'udata' have to be in the interval [0,1].")
    if (!(method %in% c("MR", "beta", "T", "TLL1", "TLL2", "TTPI", "TTCV")))
        stop("method not implemented")
    if (mult <= 0)
        stop("'mult' has to be a positive number.")
    if (is.na(knots))
        knots <- 100/d
    knots <- round(knots)
    stopifnot(is.numeric(renorm.iter))
    renorm.iter <- round(renorm.iter)
    stopifnot(renorm.iter >= 0)
    stopifnot(is.logical(info))
    
    ## bandwidth selection and adjustment (with bandwidth checks)
    if (missing(bw))
        bw <- bw_select(udata, method)
    if (any(is.na(bw)))
        bw <- bw_select(udata, method)
    if (method %in% c("TLL1", "TLL2")) {
        if (is.null(bw$B) | is.null(bw$alpha))
            stop("For methods 'TLL1/2', you have to provide a list 'bw = list(B = <your.B>,  alpha = <your.alpha>)'
where both parts of the bandwidth specification are provided via  'your.B', and 'your.alpha'.")
        if (any(c(diag(bw$B), bw$alpha) <= 0))
            stop("Bandwidths have to be positive.")
        bw$alpha <- mult * bw$alpha
    } else if (method %in% c("TTPI", "TTCV")) {
        if (bw[1] < 0)
            stop("The smoothing parameter (bw[1]) has to be positive.")
        if (bw[3] < 0)
            stop("The first tapering parameter (bw[3]) has to be positive.")
        if (abs(bw[2] > 0.9999))
            stop("The correlation parameter (bw[2]) has to lie in (-1,1).")
        bw[1] <- bw[1] * mult
    } else if (method == "T") {
        B <- as.matrix(bw)
        if (nrow(B) == 1)
            bw <- diag(as.numeric(bw), d)
        bw <- mult * bw
        bw <- bw * sqrt(min((n^(-1/6))^2 / det(bw), 1))
    } else if (method == "MR") {
        bw <- min(bw * mult, 1)
    } else if (method == "beta") {
        bw <- bw * mult
    }
    
    ## fit model for method TLL
    if (method %in% c("TLL1", "TLL2")) {
        zdata <-  qnorm(udata)
        lfit <- my_locfit(zdata,
                          bw$B,
                          bw$alpha,
                          deg = as.numeric(substr(method, 4, 4)))
    } else {
        lfit <- NULL
    }
    
    ## construct grid with k knots in dimension d
    pnts <- pnorm(seq(-3.25, 3.25, l = knots))
    grid <- as.matrix(do.call(expand.grid,
                              split(rep(pnts, d), rep(c(1, 2), each = knots))))
    
    ## evaluate estimator on grid
    evalf <- eval_func(method)  # get evaluation function
    object <- list(udata = udata,
                   bw = bw,
                   lfit = lfit,
                   method = method)
    vals <- array(evalf(grid, obj = object), dim = rep(knots, d))
    
    ## rescale copula density to have uniform margins
    if (renorm.iter > 0) {
        tmplst <- split(rep(seq.int(knots)-1, d-1),
                        ceiling(seq.int(knots*(d-1))/knots))
        helpind <- as.matrix(do.call(expand.grid, tmplst))
        vals <- renorm(vals, pnts, renorm.iter, helpind)
    }
    
    ## store results
    res <- list(udata    = udata,
                grid     = pnts,
                estimate = vals,
                bw       = bw,
                method   = method)
    class(res) <- "kdecopula"
    
    ## add effp for TLL
    if (method %in% c("TLL1", "TLL2"))
        res$effp <- eff_num_par(udata, likvalues, bw, method, lfit)
    
    ## add further information if asked for
    if (info) {
        # log-likelihood
        likvalues <- dkdecop(udata, res)
        loglik <- sum(log(likvalues))
        if (!(method %in% c("TTPI", "TTCV"))) {
            # effective number of parameters
            effp <- eff_num_par(udata, likvalues, bw, method, lfit)
            # information criteria
            AIC  <- - 2 * loglik + 2 * effp
            cAIC <- AIC + (2 * effp * (effp + 1)) / (n - effp - 1)
            BIC  <- - 2 * loglik + log(n) * effp
        } else {
#             warning("Effective number of parameters not yet implemented for this method.
# Use 'info = FALSE' if you don't want to see this message.")
            effp <- AIC <- cAIC <- BIC <- NA
        }
        
        ## store results
        res$info <- list(likvalues = likvalues,
                         loglik    = loglik,
                         effp      = effp ,
                         AIC       = AIC,
                         cAIC      = cAIC,
                         BIC       = BIC)
        
        ## remove effp for TLL
        if (method %in% c("TLL1", "TLL2"))
            res$effp <- NULL
    }
    
    ## return results as kdecopula object
    res
}
