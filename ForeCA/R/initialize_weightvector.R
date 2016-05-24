#' @title Initialize weightvector for iterative ForeCA algorithms
#' @description
#' \code{initialize_weightvector} returns a unit norm (in \eqn{\ell^2})
#' vector \eqn{\mathbf{w}_0 \in R^K} that can be used as the starting
#' point for any iterative ForeCA algorithm, e.g., 
#' \code{\link{foreca.EM.one_weightvector}}. Several 
#' quickly computable heuristics are available via the \code{method} argument. 
#' 
#' @keywords manip
#' @inheritParams common-arguments
#' @inheritParams sfa
#' @param num.series positive integer; number of time series \eqn{K} (determines the length 
#' of the weightvector). If \code{num.series = 1} it simply returns
#' a 1 \eqn{\times} 1 array equal to \code{1}.
#' @param method string; which heuristics should be used to generate a good starting \eqn{\mathbf{w}_0}?
#' Default: \code{"rnorm"}; see Details.
#' @param seed non-negative integer; seed for random initialization which will be 
#' returned for reproducibility. By default it sets a random seed.
#' @details
#' The \code{method} argument specifies the heuristics that is used to get a good
#' starting vector \eqn{\mathbf{w}_0}:
#' 
#' \itemize{
#'  \item{\code{"max"}}{ vector with all \eqn{0}s, but a \eqn{1} at the position
#'  of the maximum forecastable series in \code{U}.}
#'  \item{\code{"rcauchy"}}{ random start using \code{rcauchy(k)}.}
#'  \item{\code{"rnorm"}}{ random start using \code{rnorm(k, 0, 1)}.}
#'  \item{\code{"runif"}}{ random start using \code{runif(k, -1, 1)}.}
#'  \item{\code{"PCA.large"}}{ first eigenvector of PCA (largest variance signal).}
#'  \item{\code{"PCA.small"}}{ last eigenvector of PCA (smallest variance signal).}
#'  \item{\code{"PCA"}}{ checks both small and large, and chooses the one with higher
#'             forecastability as computed by \code{\link{Omega}}..}
#'  \item{\code{"SFA.fast"}}{ last eigenvector of SFA (fastest signal).}
#'  \item{\code{"SFA.slow"}}{ first eigenvector of SFA (slowest signal).}
#'  \item{\code{"SFA"}}{ checks both slow and fast, and chooses the one with higher
#'             forecastability as computed by \code{\link{Omega}}.}
#' } 
#' 
#' Each vector has length K and is automatically normalized to have unit norm 
#' in \eqn{\ell^2}.
#' 
#' For the \code{'SFA*'} methods see \code{\link{sfa}}. 
#' Note that maximizing (or minimizing) the lag \eqn{1} auto-correlation does 
#' not necessarily yield the most forecastable signal, but it's a good start.
#' @return 
#' numeric; a vector of length \eqn{K} with unit norm in \eqn{\ell^2}.
#' @examples
#' XX <- diff(log(EuStockMarkets))
#' \dontrun{
#' initialize_weightvector(U = XX, method = "SFA")
#' }
#' initialize_weightvector(num.series = ncol(XX), method = "rnorm")
#' @export
#' 
initialize_weightvector <- function(U = NULL, f.U = NULL, 
                                    num.series = ncol(U),
                                    method = c("rnorm", "max", "SFA", "PCA", 
                                               "rcauchy", "runif", "SFA.slow", "SFA.fast",
                                               "PCA.large", "PCA.small"), 
                                    seed = sample(1e6, 1), ...) {
        
    stopifnot(is.null(U) || is.array(U) || is.ts(U) || is.matrix(U),
              is.null(f.U) || is.array(f.U),
              num.series >= 1,
              is.character(method))
    
    if (!is.null(U)) {
      # if U is provided it must be whitened
      U <- check_whitened(U)
    }
    
    method <- match.arg(method)
    
    if (is.null(num.series) && is.null(f.U)) {
      stop("You must provide either 'num.series', 'series', or 'f.U'.")
    }
    
    if (method %in% c("SFA", "SFA.slow", "SFA.fast",
                      "PCA", "PCA.large", "PCA.small") && is.null(U)) {
      stop("For SFA- or PCA-type methods you must provide data via the 'U' argument.")
    } 
    
    if (is.null(num.series)) {
      if (is.null(f.U)) {
        num.series <- ncol(U)
      } else if (is.null(U)) {
        num.series <- dim(f.U)[2]
      } else {
        stop("Something went wrong in initialize_weightvector().")
      }
    }
    
    if (num.series == 1) {
      return(cbind(1))
    }

    set.seed(seed)
    if (method == "rcauchy") {
      ww0 <- rcauchy(num.series)
    } else if (method == "runif") {
      ww0 <- runif(num.series, -1, 1)
    } else if (method == "rnorm") {
      ww0 <- rnorm(num.series)
    } else if (method == "max") {
      ww0 <- rep(0, num.series)
      Omega.tmp <- apply(get_spectrum_from_mvspectrum(f.U), 2, 
                         function(x) Omega(mvspectrum.output = x))
      ww0[which.max(Omega.tmp)] <- 1
    } else if (any(method == c("SFA", "SFA.slow", "SFA.fast"))) {
      sfa.est <- ForeCA::sfa(U, ...)
      ww.slow <- sfa.est$loadings[, 1]
      ww.fast <- sfa.est$loadings[, num.series]
      
      if (method == "SFA.slow") {
        ww0 <- ww.slow
      } else if (method == "SFA.fast") {
        ww0 <- ww.fast
      } else if (method == "SFA") {
        stopifnot(!is.null(f.U))
        
        omega.slow <- Omega(mvspectrum.output = foreca.EM.E_step(f.U, ww.slow))
        omega.fast <- Omega(mvspectrum.output = foreca.EM.E_step(f.U, ww.fast))
        if (omega.fast > omega.slow) {
          ww0 <- ww.fast
        } else {
          ww0 <- ww.slow
        }
      }
    } else if (method %in% c("PCA", "PCA.large", "PCA.small")) {
      pca.est <- princomp(U, ...)
      ww.large <- pca.est$loadings[, 1]
      ww.small <- pca.est$loadings[, num.series]
      
      if (method == "PCA.large") {
        ww0 <- ww.large
      } else if (method == "PCA.small") {
        ww0 <- ww.small
      } else if (method == "PCA") {
        stopifnot(!is.null(f.U))
        
        omega.large <- Omega(mvspectrum.output = foreca.EM.E_step(f.U, ww.large))
        omega.small <- Omega(mvspectrum.output = foreca.EM.E_step(f.U, ww.small))
        if (omega.small > omega.large) {
          ww0 <- ww.small
        } else {
          ww0 <- ww.large
        }
      }
    } else {
      stop("Method '", method, "' is not available for initialize_weightvector().")
    }
    
    if (ww0[1] != 0) {
      # make the first entry always positive (if it's not zero)
      ww0 <- ww0 * sign(ww0[1])
    }
    # normalize
    ww0 <- rbind(ww0 / base::norm(ww0, "2"))
    rownames(ww0) <- NULL
    return(ww0)
  } 
