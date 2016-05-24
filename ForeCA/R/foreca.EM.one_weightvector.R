#' @title EM-like algorithm to estimate optimal ForeCA transformation
#' @name foreca.EM.one_weightvector
#' @keywords manip optimize iteration
#' @description
#' \code{foreca.EM.one_weightvector} finds the optimal weightvector \eqn{\mathbf{w}^*} 
#' that gives the most forecastable signal \eqn{y_t^* = \mathbf{U}_t \mathbf{w}^*} 
#' using an EM-like algorithm (see References).
#' @inheritParams common-arguments
#' @param ... other arguments passed to \code{\link{mvspectrum}}
#' @param init.weightvector numeric; starting point \eqn{\mathbf{w}_0} for several
#' iterative algorithms.  By default it uses a (normalized) random vector from a 
#' standard Normal distribution (see \code{\link{initialize_weightvector}}).
#' @export
#' @references 
#' Goerg, G. M. (2013). \dQuote{Forecastable Component Analysis}. 
#' Journal of Machine Learning Research (JMLR) W&CP 28 (2): 64-72, 2013.
#' Available at \url{jmlr.org/proceedings/papers/v28/goerg13.html}.
#' @seealso
#' \code{\link{foreca.one_weightvector}}, \code{\link{foreca.EM-aux}}
#' @return
#' A list with useful quantities like the optimal weighvector, the corresponding
#' signal, and its forecastability. 
#' @examples
#' \dontrun{
#' XX <- diff(log(EuStockMarkets)[100:200,]) * 100
#' one.weight <- foreca.EM.one_weightvector(whiten(XX)$U, 
#'                                          spectrum.control = 
#'                                             list(method = "wosa"))
#' }
#' 
foreca.EM.one_weightvector <- function(U, f.U = NULL, 
                                       spectrum.control = list(),
                                       entropy.control = list(),
                                       algorithm.control = list(),
                                       init.weightvector = 
                                         initialize_weightvector(num.series = ncol(U),
                                                                 method = 'rnorm'),
                                       ...) {

  UU <- U
  if (!is.ts(UU)) {
    UU <- ts(UU)
  }
  
  num.series <- ncol(UU)
  num.freqs <- floor(nrow(UU) / 2)
  
  spectrum.control <- complete_spectrum_control(spectrum.control)
  algorithm.control <- complete_algorithm_control(algorithm.control)
  
  entropy.control$base <- NULL
  entropy.control <- complete_entropy_control(entropy.control, 
                                              num.outcomes = 2 * num.freqs)
  
  # check that init.weightvector has same dimension as series
  stopifnot(length(init.weightvector) == num.series)
  if (!isTRUE(all.equal(target = 1,
                        current = base::norm(init.weightvector, "2")))) {
    warning("Initial weightvector did not have unit norm. It was normalized automatically.")
    init.weightvector <- init.weightvector / base::norm(init.weightvector, "2")
  }

  if (is.null(f.U)) {
    f.U <- mvspectrum(UU, method = spectrum.control$method, normalize = TRUE, 
                        ...)
    if (!is.null(spectrum.control$kernel)) {
      f.U <- sweep(f.U, 1, spectrum.control$kernel(attr(f.U, "frequency")), 
                   FUN = "*")
      f.U <- normalize_mvspectrum(f.U)
    }
  }
  converged <- FALSE  
  wv.trace <- matrix(NA, ncol = num.series, 
                     nrow = algorithm.control$max.iter + 1)  # add one for last recording
  wv.trace[1, ] <- init.weightvector
  
  spec.ent.trace <- rep(NA, algorithm.control$max.iter + 1)
  spec.ent.trace.univ <- spec.ent.trace 
  Omega.trace <- spec.ent.trace
  Omega.trace.univ <- Omega.trace
  
  warning.msg <- NULL

  for (iter in seq_len(algorithm.control$max.iter)) {
    ww.current <- wv.trace[iter, ]
    yy.current <- c(UU %*% ww.current)
    if (spectrum.control$smoothing) {
      f.current <- mvspectrum(yy.current,
                              method = spectrum.control$method, 
                              normalize = TRUE, smoothing = TRUE)
    } else {
      f.current <- foreca.EM.E_step(f.U = f.U,
                                    weightvector = ww.current)
    }
    spec.ent.trace.univ[iter] <- 
      spectral_entropy(mvspectrum.output = f.current,
                       entropy.control = entropy.control)
    spec.ent.trace[iter] <- 
      foreca.EM.h(weightvector.new = wv.trace[iter, ], 
                  f.U = f.U, f.current = f.current, 
                  entropy.control = entropy.control)
    
    Omega.trace[iter] <- Omega(mvspectrum.output = f.current,
                               entropy.control = entropy.control)
    Omega.trace.univ[iter] <- Omega(series = yy.current,
                                      entropy.control = entropy.control,
                                      spectrum.control = spectrum.control)
    if (converged) {
      break
    }
    wv.trace[iter + 1,] <- foreca.EM.M_step(f.U, f.current, 
                                            minimize = TRUE,
                                            entropy.control = 
                                              entropy.control)$vector
    if (iter > 1) {
      # round to avoid numerical rounding errors
      abs.change <- spec.ent.trace[iter] - spec.ent.trace[iter - 1]  # this is <0
      if (abs.change > 0) {
        warning.msg <- paste0("Spectral entropy increased (!) in iteration ",
                              iter, " of foreca.EM.one_weightvector.\n ",
                              "If this is not the last iteration, please check results.")
        warning.iter <- iter
      } else {
        # decrease from previous iteration
        # rel.change <- spec.ent.trace[iter] / spec.ent.trace[iter - 1] - 1
        if (abs(abs.change) < algorithm.control$tol) {
          # convergence stop
          converged <- TRUE  # break in next iteration
        }
      }
    }
  }  # end of iter
  
  if (iter == algorithm.control$max.iter) {
    spec.ent.trace[iter + 1] <- spec.ent.trace[iter]
  }
  wv.trace <- na.omit(wv.trace)
  spec.ent.trace <- na.omit(spec.ent.trace)
  spec.ent.trace.univ <- na.omit(spec.ent.trace.univ)
  Omega.trace.univ <- na.omit(Omega.trace.univ)
  Omega.trace <- na.omit(Omega.trace)
  attr(Omega.trace, "na.action") <- NULL
  attr(Omega.trace.univ, "na.action") <- NULL
  attr(spec.ent.trace.univ, "na.action") <- NULL
  attr(wv.trace, "na.action") <- NULL
  attr(spec.ent.trace, "na.action") <- NULL
  
  wv.trace <- sweep(wv.trace, 1, sign(wv.trace[, 1]), "*")
  
  rownames(wv.trace) <- paste("Iter", seq_len(nrow(wv.trace)) - 1)
  colnames(wv.trace) <- colnames(UU)
  
  wv.final <- tail(wv.trace, 1)
  rownames(wv.final) <- NULL
  
  out <- list(weightvector.trace = wv.trace,
              h.trace = spec.ent.trace,
              h.trace.univ = spec.ent.trace.univ,
              Omega.trace = Omega.trace,
              Omega.trace.univ = Omega.trace.univ,
              Omega = tail(Omega.trace, 1),
              h = tail(spec.ent.trace, 1),
              weightvector = wv.final,
              iterations = nrow(wv.trace) - 1,
              converged = converged)
  
  if (!is.null(warning.msg) && warning.iter < out$iterations) {
    warning(warning.msg)
  }
  
  class(out) <- "foreca.EM.one_weightvector"
  invisible(out)
}