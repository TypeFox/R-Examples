#' @title Forecastable Component Analysis
#' @name foreca
#' @description 
#' \code{foreca} performs Forecastable Component Analysis (ForeCA) on 
#' \eqn{\mathbf{X}_t} -- a \eqn{K}-dimensional time series with \eqn{T} 
#' observations. Users should only call
#' \code{foreca}, rather than \code{foreca.one_weightvector} or 
#' \code{foreca.multiple_weightvectors}.
#' 
#' @inheritParams common-arguments
#' @param n.comp positive integer; number of components to be extracted. 
#' Default: \code{2}.
#' @param ... additional arguments passed to available ForeCA algorithms.
#' 
#' @return
#' An object of class \code{foreca}, which is similar to the output from \code{\link[stats]{princomp}}, 
#' with the following components (amongst others):
#' \itemize{
#' \item{\code{center}:}{ sample mean \eqn{\widehat{\mu}_X} of each \code{series},}
#' \item{\code{whitening}:}{ whitening matrix of size \eqn{K \times K} 
#' from \code{\link{whiten}}: \eqn{\mathbf{U}_t = (\mathbf{X}_t - \widehat{\mu}_X) \cdot whitening};
#' note that \eqn{\mathbf{X}_t} is centered prior to the whitening transformation,}
#' \item{\code{weightvectors}:}{ orthonormal matrix of size \eqn{K \times n.comp},
#' which converts whitened data to \code{n.comp} forecastable components (ForeCs) 
#' \eqn{\mathbf{F}_t = \mathbf{U}_t \cdot weightvectors}, }
#' \item{\code{loadings}:}{ combination of whitening \eqn{\times} weightvectors to obtain the final
#' loadings for the original data: 
#' \eqn{\mathbf{F}_t = (\mathbf{X}_t - \widehat{\mu}_X) \cdot whitening \cdot 
#' weightvectors}; again, it centers \eqn{\mathbf{X}_t} first,}
#' \item{\code{loadings.normalized}:}{ normalized loadings (unit norm).  Note
#' though that if you use these normalized loadings the resulting
#' signals do not have variance 1 anymore.}
#' \item{\code{scores}:}{ \code{n.comp} forecastable components \eqn{\mathbf{F}_t}. 
#'  They have mean 0, variance 1, and are uncorrelated.}
#' \item{\code{Omega}:}{ forecastability score of each ForeC of \eqn{\mathbf{F}_t}.}
#' }
#' 
#' ForeCs are ordered from most to least forecastable (according to 
#' \code{\link{Omega}}).
#' 
#' @section Warning:
#' 
#' Estimating Omega directly from the ForeCs \eqn{\mathbf{F}_t} can be different 
#' to the reported \code{$Omega} estimates from \code{foreca}.  Here is why:
#' 
#' In theory \eqn{f_y(\lambda)} of a linear combination 
#' \eqn{y_t = \mathbf{X}_t \mathbf{w}} can be analytically computed from 
#' the multivariate spectrum \eqn{f_{\mathbf{X}}(\lambda)} by the 
#' quadratic form
#' \eqn{f_y(\lambda) = \mathbf{w}' f_{\mathbf{X}}(\lambda) \mathbf{w}} for all 
#' \eqn{\lambda} (see \code{\link{spectrum_of_linear_combination}}).
#' 
#' In practice, however, this identity does not hold always exactly since 
#' (often data-driven) control setting for spectrum estimation are not identical
#' for the high-dimensional, noisy
#' \eqn{\mathbf{X}_t} and the combined univariate time series \eqn{y_t} 
#' (which is usually more smooth, less variable). Thus estimating
#' \eqn{\widehat{f}_y} directly  from \eqn{y_t} can give slightly different
#' estimates to computing it as \eqn{\mathbf{w}'\widehat{f}_{\mathbf{X}}\mathbf{w}}.  Consequently also \code{Omega} estimates
#' can be different. 
#' 
#' 
#' In general, these differences are small and have no relevant implications
#' for estimating ForeCs.  However, especially for rare occasions, the obtained ForeCs can have
#' smaller \code{Omega} than the maximum \code{Omega} of the original series.
#' In such a case users should not re-estimate \eqn{\Omega} from the resulting
#' ForeCs \eqn{\mathbf{F}_t}, but access them via \code{$Omega} provided 
#' by \code{'foreca'} output (the univariate estimates are stored in \code{$Omega.univ}).
#' 
#' @references
#' Goerg, G. M. (2013). \dQuote{Forecastable Component Analysis}. 
#' Journal of Machine Learning Research (JMLR) W&CP 28 (2): 64-72, 2013.
#' Available at \url{jmlr.org/proceedings/papers/v28/goerg13.html}.
#' @export
#' @examples
#' XX <- diff(log(EuStockMarkets[c(100:200),])) * 100
#' plot(ts(XX))
#' \dontrun{
#' ff <- foreca(XX[,1:4], n.comp = 2, plot = TRUE)
#' ff
#' summary(ff)
#' plot(ff)
#' }
#' 

foreca <- function(series, n.comp = 2, algorithm.control = list(type = "EM"),
                   ...) {
  
  algorithm.control <- complete_algorithm_control(algorithm.control)
  series.name <- deparse(substitute(series))
  
  num.series <- ncol(series)
  if (is.null(num.series) || num.series == 1) {
    stop("You must provide at least 2 time series (columns).")
  }
  # number of components must be smaller or equal to number of time series
  if (n.comp > num.series) {
    stop("You can not extract ", n.comp, " components from only ", num.series, " time series.\n",
         "Please reduce n.comp to at most ", num.series, ".")
  }
  
  PW.all <- whiten(series)
  UU <- PW.all$U
  
  if (algorithm.control$type == "EM"){
    out <- foreca.multiple_weightvectors(UU, 
                                         algorithm.control = algorithm.control, 
                                         n.comp = n.comp, 
                                         dewhitening = PW.all$dewhitening, 
                                         ...)
  } else {
    stop("Algorithm type '", algorithm.control$type, "' is not implemented.")
  }
  out$weightvectors <- cbind(out$weightvectors)
  out$whitening <- PW.all$whitening
  out$loadings <- PW.all$whitening %*% out$weightvectors
  
  rownames(out$loadings) <- colnames(series)
  colnames(out$loadings) <- paste0("ForeC", seq_len(ncol(out$loadings)))
  class(out$loadings) <- "loadings"

  out <- c(out,
           list(n.obs = nrow(series),
                Omega.matrix = 
                  out$loadings %*% diag(out$Omega, 
                                        nrow = length(out$Omega)) %*% 
                  t(out$loadings),
                Omega.U.matrix = out$weightvectors %*% diag(out$Omega,
                                                            nrow = length(out$Omega)) %*% 
                  t(out$weightvectors),
                sdev = out$Omega,
                lambdas = out$Omega,
                center = PW.all$center,
                # remove mean from each series so they become zero mean
                # series.centered = sweep(series, 2, colMeans(series), "-"),
                series = series))
  out$loadings.normalized <- sweep(out$loadings, 2, 
                                   base::norm(out$loadings, "2"), FUN = "/")
  out$series.name <- series.name
  class(out) <- c("foreca", class(out))
  return(out)
}

#' @rdname foreca
#' @description
#' \code{foreca.one_weightvector} is a wrapper around several algorithms that 
#' solve the ForeCA optimization problem for a single weightvector \eqn{\mathbf{w}_i}
#' and whitened time series \eqn{\mathbf{U}_t}.
#' @param keep.all.optima logical; if \code{TRUE}, it keeps the optimal 
#' solutions of each random start. Default: \code{FALSE} (only returns the best solution).
#' @param dewhitening optional; if provided (returned by \code{\link{whiten}})
#' then it uses the dewhitening transformation to obtain the original 
#' series \eqn{\mathbf{X}_t} and it uses that vector (normalized) as the initial 
#' weightvector which corresponds to the series \eqn{\mathbf{X}_{t,i}} 
#' with larges \code{\link{Omega}}.
#' @examples
#'
#' PW <- whiten(XX)
#' one.weight.em <- foreca.one_weightvector(U = PW$U,
#'                                         dewhitening = PW$dewhitening,
#'                                         algorithm.control = 
#'                                           list(num.starts = 2,
#'                                                type = "EM"),
#'                                         spectrum.control = 
#'                                           list(method = 'wosa'))
#' plot(one.weight.em)
#' 
#' @export

foreca.one_weightvector <- function(U, f.U = NULL, 
                                    spectrum.control = list(),
                                    entropy.control = list(),
                                    algorithm.control = list(),
                                    keep.all.optima = FALSE,
                                    dewhitening = NULL,
                                    ...) {
  
  UU <- check_whitened(U)
  num.series <- ncol(UU)
  num.obs <- nrow(UU)
  
  stopifnot(is.matrix(UU) || is.ts(UU),
            num.obs > 1,
            num.series > 1,
            class(f.U) == "mvspectrum" || is.null(f.U),
            is.logical(keep.all.optima))
  
  if (!is.null(dewhitening)) {
    stopifnot(is.matrix(dewhitening),
              num.series == nrow(dewhitening))
  }
  
  if (!is.ts(UU)) {
    UU <- ts(UU)
  }
  
  spectrum.control <- complete_spectrum_control(spectrum.control)
  algorithm.control <- complete_algorithm_control(algorithm.control)
  
  num.freqs <- floor(num.obs / 2)
  entropy.control$base <- NULL
  entropy.control <- complete_entropy_control(entropy.control, 
                                              num.outcomes = 2 * num.freqs)
  
  if (is.null(f.U)) {
    f.U <- mvspectrum(UU, method = spectrum.control$method, normalize = TRUE, 
                      ...)
    if (!is.null(spectrum.control$kernel)) {
      f.U <- sweep(f.U, 1, spectrum.control$kernel(attr(f.U, "frequency")), 
                   FUN = "*")
      f.U <- normalize_mvspectrum(f.U)
    }
  }
  
  check_mvspectrum_normalized(f.U)
  
  Omega.best <- 0
  converged <- FALSE
  
  # Prepare arguments for any one vector algorithm
  args.for.one_weightvector.algo <- 
    list(spectrum.control = spectrum.control,
         entropy.control = entropy.control,
         algorithm.control = algorithm.control,
         U = UU,
         f.U = f.U)
  
  all.optima <- list()
  
  init.methods <- c("rnorm", "SFA.slow", "max", "SFA.fast", "rcauchy", "runif")
  num.init.methods <- length(init.methods)
  
  init.weights.matrix <- matrix(NA, ncol = num.series, nrow = 1)
  method.names <- c()
  
  for (TRY in seq_len(algorithm.control$num.starts)) {
    if (TRY <= num.init.methods) {
      current.method <- init.methods[TRY]
      init.weights <- initialize_weightvector(UU, f.U, method = current.method)
    } else if (TRY == num.init.methods + 1) {
      current.method <- paste("SFA with lag", frequency(UU), "difference.")
      init.weights <- initialize_weightvector(UU, f.U, method = "SFA.slow", 
                                              lag = frequency(UU))
    } else if (TRY == num.init.methods + 2) {
      current.method <- "Whitened series with highest Omega."
      Omega.U <- Omega(UU, spectrum.control = spectrum.control,
                       entropy.control = entropy.control)
      init.weights <- rep(0, num.series)
      init.weights[which.max(Omega.U)] <- 1
    } else if (TRY == num.init.methods + 2) {
      current.method <- "Original series with highest Omega."
    } else { 
      current.method <- "rnorm"
      # do random draw from normal distribution every run
      init.weights <- initialize_weightvector(UU, f.U, method = "rnorm")
    }
    init.weights.matrix <- rbind(init.weights.matrix, 
                                 init.weights)
    method.names <- c(method.names, current.method) 
  }
  # remove first NA row
  init.weights.matrix <- init.weights.matrix[-1,]
  
  if (!is.null(dewhitening)) {
    orig.series <- UU %*% dewhitening
    Omega.orig <- Omega(orig.series, 
                        spectrum.control = spectrum.control,
                        entropy.control = entropy.control)
    init.weights <- dewhitening[, which.max(Omega.orig)]
    init.weights <- init.weights / base::norm(init.weights, "2")
    current.method <- "Original series with highest Omega."
    
    init.weights.matrix <- rbind(init.weights.matrix, 
                                 init.weights)
    method.names <- c(method.names, current.method) 
  }
  
  if (init.weights[1] != 0) {
    # make the sign the same for first entry
    init.weights.matrix <- t(apply(init.weights.matrix, 1,
                                   function(x) {
                                     if (x[1] != 0) {
                                       return(x / sign(x)[1])
                                     } else {
                                       return(x)
                                     }}))
  }
  
  for (ii in seq_len(nrow(init.weights.matrix))) {
    # pass init.weights to control settings for algorithm
    args.for.one_weightvector.algo$init.weightvector <- 
      init.weights.matrix[ii,]
    
    init.series <- UU %*% args.for.one_weightvector.algo$init.weightvector
    Omega.init.tmp <- Omega(init.series, 
                            spectrum.control = spectrum.control,
                            entropy.control = entropy.control)
    if (algorithm.control$type == "EM") {
      one.results.tmp <- do.call(foreca.EM.one_weightvector, 
                                 args.for.one_weightvector.algo)
    } else {
      stop("Algorithm '", algorithm.control$type, "' is not implemented.")
    }
    Omega.tmp <- one.results.tmp$Omega
    if (keep.all.optima) {
      all.optima[[ii]] <- one.results.tmp
    }
    
    if (Omega.tmp > Omega.best) {
      best.attempt <- method.names[ii]
      one.results.best <- one.results.tmp
      Omega.init <- Omega.init.tmp
      Omega.best <- Omega.tmp
      converged <- one.results.tmp$converged
    }
  }  # end of ii
  
  if (!converged) {
    warning("Convergence has not been reached. Please try again.")
  }
  
  wv.trace.best <- one.results.best$weightvector.trace
  rownames(wv.trace.best) <- paste("Iter", seq_len(nrow(wv.trace.best)) - 1)
  colnames(wv.trace.best) <- colnames(UU)
  
  one.results.best$weightvector.trace <- 
    .make_smooth_weightvector_trace(one.results.best$weightvector.trace)
  
  out <- list(estimate = one.results.best,
              algorithm.control = algorithm.control,
              spectrum.control = spectrum.control,
              entropy.control = entropy.control,
              best.attempt = best.attempt)
  
  out <- c(out,
           list(h = out$estimate$h,
                iterations = out$estimate$iterations,
                Omega.trace = out$estimate$Omega.trace,
                converged = converged,
                weightvector = 
                  as.matrix(tail(out$estimate$weightvector.trace, 1))))
  rownames(out$weightvector) <- NULL
  
  out <- c(out,
           list(score = ts(UU %*% t(out$weightvector)),
                Omega.init = Omega.init,
                Omega = Omega.best,
                best.f = 
                  spectrum_of_linear_combination(f.U, out$weightvector)))
  out$score <- check_whitened(out$score, FALSE)
  out$best.f.univ <- c(mvspectrum(out$score,
                                  spectrum.control = spectrum.control,
                                  normalize = TRUE))
  out$Omega.univ <- Omega(out$score,
                          spectrum.control = spectrum.control,
                          entropy.control = entropy.control)
  if (keep.all.optima) {
    out$all.optima <- all.optima
  }
  class(out) <- "foreca.one_weightvector"
  invisible(out)
} 

# This function multiplies the rows of the trace matrix by +/- 1 so 
# that the path becomes a smooth curve (and it not jumping back and forth).
.make_smooth_weightvector_trace <- function(weightvector.trace) {
   
  num.iter <- nrow(weightvector.trace)
  # TODO: do this in closed form given the smoothness
  for (ii in seq_len(num.iter)) {
    iter.smoothness <- c(0, apply(diff(weightvector.trace), 1, 
                                  base::norm, "2"))
    if (any(iter.smoothness > 1)) {
      first.switch <- which(iter.smoothness > 1)[1]
      multipliers <- rep(1, length = num.iter)
      multipliers[first.switch] <- -1
      weightvector.trace <- sweep(weightvector.trace, 1, multipliers, "*")
    } else {
      # if all differences are smooth; then break and return
      break
    }
  }
  return(weightvector.trace)
}


#' @rdname foreca
#' @description
#' \code{foreca.multiple_weightvectors} applies \code{foreca.one_weightvector}
#' iteratively to \eqn{\mathbf{U}_t} in order to obtain multiple weightvectors
#' that yield most forecastable, uncorrelated signals.
#' @inheritParams foreca
#' @param plot logical; if \code{TRUE} a plot of the current optimal 
#' solution \eqn{\mathbf{w}_i^*} will be shown and updated for each iteration 
#' \eqn{i = 1, ..., } \code{n.comp} of any iterative algorithm. Default: \code{FALSE}.
#' @export
#' @keywords iteration
#' @examples
#' \dontrun{
#' 
#' PW <- whiten(XX)
#' ff <- foreca.multiple_weightvectors(PW$U, n.comp = 2,
#'                                     dewhitening = PW$dewhitening)
#' ff
#' plot(ff$scores)
#' }

foreca.multiple_weightvectors <- function(U,
                                          spectrum.control = list(),
                                          entropy.control = list(),
                                          algorithm.control = list(),
                                          n.comp = 2,
                                          plot = FALSE,
                                          dewhitening = NULL,
                                          ...) {
  
  UU <- U
  UU <- check_whitened(UU)
  
  num.series <- ncol(UU)
  num.obs <- nrow(UU)
  # if multiple methods are specified (default), it uses the first [1] by default
  spectrum.control <- complete_spectrum_control(spectrum.control)
  
  num.freqs <- floor(num.obs / 2)
  entropy.control$base <- NULL
  entropy.control <- complete_entropy_control(entropy.control,
                                              num.outcomes = 2 * num.freqs)
  
  algorithm.control <- complete_algorithm_control(algorithm.control)
  
  out <- list(spectrum.control = spectrum.control,
              entropy.control = entropy.control,
              algorithm.control = algorithm.control)
  
  out <- c(out,
           list(weightvectors = matrix(NA, ncol = n.comp, nrow = num.series),
                scores = matrix(NA, ncol = n.comp, nrow = num.obs),
                h = rep(NA, n.comp),
                Omega = rep(NA, n.comp),
                Omega.univ = rep(NA, n.comp),
                weights = list()))
  
  null.space <- diag(1, num.series)
  for (comp.id in seq_len(n.comp)) {
    UU.in.null <- as.matrix(UU) %*% null.space
    attr(UU.in.null, "whitened") <- TRUE
    if (is.null(dewhitening)) {
      dewhitening.tmp <- NULL
    } else {
      dewhitening.tmp <- t(null.space) %*% dewhitening
    }
    if (comp.id == num.series) {
      out$h[comp.id] <- spectral_entropy(UU.in.null, 
                                         spectrum.control = spectrum.control, 
                                         entropy.control = entropy.control)
      out$scores[, n.comp] <- UU.in.null
      out$Omega[comp.id] <- Omega(UU.in.null,
                                  spectrum.control = spectrum.control,
                                  entropy.control = entropy.control)
      out$Omega.univ[comp.id] <- Omega(UU.in.null,
                                       spectrum.control = spectrum.control,
                                       entropy.control = entropy.control)
      out$weights[[comp.id]] <- 1
    } else {
      # This is the actual algorithm that obtains the new w*
      opt.wv <- foreca.one_weightvector(U = UU.in.null,
                                        f.U = NULL,
                                        spectrum.control = spectrum.control,
                                        entropy.control = entropy.control,
                                        algorithm.control = algorithm.control,
                                        dewhitening = dewhitening.tmp,
                                        ...)
      # plot results
      if (plot) {
        plot(opt.wv, main = paste("Component", comp.id))
      }
      out$h[comp.id] <- opt.wv$h
      out$scores[, comp.id] <- opt.wv$score
      out$weights[[comp.id]] <- opt.wv$weightvector
      out$Omega[comp.id] <- opt.wv$Omega
      out$Omega.univ[comp.id] <- opt.wv$Omega.univ
    }
    out$weightvectors[, comp.id] <- null.space %*% matrix(out$weights[[comp.id]],
                                                          ncol = 1)
    null.space <- MASS::Null(out$weightvectors[, seq_len(comp.id)])
  }
  
  omega.order <- order(out$Omega.univ, decreasing = TRUE)
  if (!identical(sort(omega.order), omega.order)) {
    warning("'foreca.multiple_weightvectors()' did not extract ForeCs in order of ",
            "decreasing forecastability; ",
            "loadings and sources were re-ordered accordingly. ",
            "Please check results.")
    cat("Original order:\n")
    cat("\t", paste(omega.order, collapse = ", "), "\n")
    
    out$scores <- out$scores[, omega.order]
    out$weightvectors <- out$weightvectors[, omega.order]
    out$Omega <- out$Omega[omega.order]
    out$Omega.univ <- out$Omega.univ[omega.order]
    out$h <- out$h[omega.order]
  }
  colnames(out$scores) <- paste0("ForeC", seq_len(ncol(out$scores)))
  out$scores <- ts(out$scores)
  class(out$weightvectors) <- "loadings"
  names(out$Omega) <- colnames(out$scores)
  
  class(out) <- "foreca.multiple_weightvectors"
  return(out)
}
