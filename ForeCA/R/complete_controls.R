#' @title Completes several control settings
#' @name complete-controls
#' @description
#' Completes algorithm, entropy, and spectrum control lists.
#' 
#' @param spectrum.control list; control settings for spectrum estimation.
#' @param entropy.control list; control settings for entropy estimation.
#' @param algorithm.control list; control parameters for any \emph{iterative} ForeCA 
#' algorithm.
#' @param num.outcomes positive integer; number of outcomes for the discrete probability
#' distribution. Must be specified (no default value).
#' @return 
#' A list with fully specified algorithm, entropy, or spectrum controls. 
#' Default values are only added if the input \code{{spectrum,entropy,algorithm}.control} 
#' list does not already set this value.
#' 
#' @seealso \code{\link{mvspectrum}}, \code{\link{discrete_entropy}}, 
#' \code{\link{continuous_entropy}}
#' @keywords utils
#' 
NULL

#' @rdname complete-controls
#' @return
#' \code{complete_algorithm_control} returns a list containing:
#' 
#' \item{max.iter}{maximum number of iterations; default: \code{50}.}
#' \item{num.starts}{number of random starts to avoid local optima; default: \code{10}.}
#' \item{tol}{tolerance for when convergence is reached in any \emph{iterative} 
#' ForeCA algorithm;  default: \code{1e-03}.}
#' \item{type}{string; type of algorithm. Default: \code{'EM'}.}
#' @export
#' 
complete_algorithm_control <- function(algorithm.control = 
                                         list(max.iter = 50,
                                              num.starts = 10,
                                              tol = 1e-3,
                                              type = 'EM')) {
  stopifnot(inherits(algorithm.control, "list"))
  
  valid.entries <- c("max.iter", "num.starts", "tol", "type")
  if (!is.null(names(algorithm.control))) {
    matched <- match(names(algorithm.control), valid.entries)
    if (any(is.na(matched))) {
      print(names(algorithm.control[is.na(matched)]))
      stop("'algorithm.control' has unvalid entries (see above). Please remove them from the list.")
    }
  }
            
  if (is.null(algorithm.control$tol)) {
    algorithm.control$tol <- 1e-6
  } 
  stopifnot(is.numeric(algorithm.control$tol),
            length(algorithm.control$tol) == 1,
            algorithm.control$tol > .Machine$double.eps^0.9)
  
  # Currently only the 'EM' type algorithm is available
  if (is.null(algorithm.control$type)) {
    algorithm.control$type <- "EM"
  }
  stopifnot(is.character(algorithm.control$type),
            length(algorithm.control$type) == 1)
  
  if (is.null(algorithm.control$max.iter)) {
    algorithm.control$max.iter <- 50
  }
  stopifnot(is.numeric(algorithm.control$max.iter),
            length(algorithm.control$max.iter) == 1,
            algorithm.control$max.iter > 0)
  
  if (is.null(algorithm.control$num.starts)) {
    algorithm.control$num.starts <- 10
  }
  stopifnot(is.numeric(algorithm.control$num.starts),
            length(algorithm.control$num.starts) == 1,
            algorithm.control$num.starts >= 0)
  return(algorithm.control)
}



#' @rdname complete-controls
#' @return
#' \code{complete_entropy_control} returns a list with:
#' 
#' \item{base}{logarithm base for the entropy.}
#' \item{method}{string; method to estimate entropy; default: \code{"MLE"}.}
#' \item{prior.probs}{prior distribution; default: uniform 
#' \code{rep(1 / num.outcomes, num.outcomes)}.}
#' \item{prior.weight}{weight of the prior distribution; default: \code{1e-3}.}
#' \item{threshold}{non-negative float; set probabilities below threshold to 
#' zero;  default: \code{0}.}
#' @export
#' 
complete_entropy_control <- function(entropy.control = 
                                       list(base = NULL, 
                                            method = "MLE",
                                            prior.probs = NULL,
                                            prior.weight = 1e-3,
                                            threshold = 0),
                                     num.outcomes) {
  stopifnot(inherits(entropy.control, "list"),
            num.outcomes > 0,
            round(num.outcomes) == num.outcomes)  # integer
  
  valid.entries <- c("base", "method", "prior.probs", "prior.weight", "threshold")
  if (!is.null(names(entropy.control))) {
    matched <- match(names(entropy.control), valid.entries)
    if (any(is.na(matched))) {
      print(names(entropy.control[is.na(matched)]))
      stop("'entropy.control' has unvalid entries (see above). Please remove them from the list.")
    }
  }
  
  # in alphabetical order
  if (is.null(entropy.control$base)) {
    entropy.control$base <- num.outcomes
  } 
  stopifnot(is.numeric(entropy.control$base),
            length(entropy.control$base) == 1,
            entropy.control$base > 0)
  
  if (is.null(entropy.control$method)) {
    entropy.control$method <- "MLE"
  }
  stopifnot(is.character(entropy.control$method),
            length(entropy.control$method) == 1)
  
  if (is.null(entropy.control$prior.probs)) {
    entropy.control$prior.probs <- rep(1 / num.outcomes, 
                                       length = num.outcomes)
  } 
  
  stopifnot(isTRUE(all.equal(1, sum(entropy.control$prior.probs))),
            all(entropy.control$prior.probs >= 0))
  #length(entropy.control$prior.probs) == num.outcomes)
  
  if (is.null(entropy.control$prior.weight)) {
    entropy.control$prior.weight <- 1e-3
  }
  stopifnot(is.numeric(entropy.control$prior.weight),
            length(entropy.control$prior.weight) == 1,
            entropy.control$prior.weight >= 0 && entropy.control$prior.weight <= 1)
  
  
  if (is.null(entropy.control$threshold)) {
    entropy.control$threshold <- 0
  }
  stopifnot(entropy.control$threshold >= 0,
            length(entropy.control$threshold) == 1)

  return(entropy.control)
}

#' @rdname complete-controls
#' @return
#' \code{complete_spectrum_control} returns a list containing:
#' 
#' \item{kernel}{R function; function to weigh each Fourier frequency \eqn{\lambda}; 
#' default: \code{NULL} (no re-weighting).}
#' \item{method}{string; method to estimate the spectrum; default: 
#' \code{'wosa'} if \pkg{sapa} is installed, \code{'mvspec'} 
#' if only \pkg{astsa} is installed, and \code{'pgram'} if
#' neither is installed.}
#' \item{smoothing}{logical; default: \code{FALSE}.}
#' 
#' Available methods for spectrum estimation are (alphabetical order)
#' 
#' \item{"ar"}{ autoregressive spectrum fit via \code{\link[stats]{spec.ar}}; 
#' only for univariate time series.}
#' \item{"direct"}{ raw periodogram using \code{\link[sapa]{SDF}}.}
#  \item{"lag window"}{ average over a neighborhood of lags in the periodogram 
#  estimator using \code{\link[sapa]{SDF}}.}
#' \item{"multitaper"}{ tapering the periodogram using \code{\link[sapa]{SDF}}.}
#' \item{"mvspec"}{ smoothed estimate using \code{\link[astsa]{mvspec}}; many tuning parameters
#' are available -- they can be passed as additional arguments (\code{...}) 
#' to \code{mvspectrum}.}
#' \item{"pgram"}{ uses \code{\link{mvpgram}}; is the same as the 
#' \code{'direct'} method, but does not rely on the \code{\link[sapa]{SDF}} 
#' package.}
#' \item{"wosa"}{ Welch overlapping segment averaging (WOSA) using \code{\link[sapa]{SDF}}.}
#' 
#' Setting \code{smoothing  = TRUE} will smooth the estimated spectrum
#' (again); this option is only available for univariate time series/spectra.
#' 
#' @export

complete_spectrum_control <- function(spectrum.control = 
                                        list(kernel = NULL, 
                                             method = c("wosa", "direct", "multitaper", 
                                                        "mvspec", "ar", "pgram"),
                                             smoothing = FALSE)) {
  stopifnot(inherits(spectrum.control, "list"))
  
  valid.entries <- c("kernel", "method", "smoothing", "taper")
  if (!is.null(names(spectrum.control))) {
    matched <- match(names(spectrum.control), valid.entries)
    if (any(is.na(matched))) {
      print(names(spectrum.control[is.na(matched)]))
      stop("'spectrum.control' has unvalid entries (see above). Please remove them from the list.")
    }
  }
  
  if (is.null(spectrum.control$kernel)) {
    spectrum.control$kernel <- NULL
  } else {
    if (!is.function(spectrum.control$kernel)) {
      stop("'kernel' must be an R function or NULL.")
    }
  }
  
  if (is.null(spectrum.control$method)) {
    if (requireNamespace("sapa", quietly = TRUE)) {
      spectrum.control$method <- "wosa"
    } else if (requireNamespace("astsa", quietly = TRUE)) {
      spectrum.control$method <- "mvspec"
    } else {
      spectrum.control$method <- "pgram"
    }
  } else if (length(spectrum.control$method) > 1) {
    # take the first method if more than one is specified
    spectrum.control$method <- spectrum.control$method[1]
  }
  stopifnot(is.character(spectrum.control$method),
            length(spectrum.control$method) == 1)
  if (spectrum.control$method %in% c("wosa", "multitaper", "direct")) {
    if (!requireNamespace("sapa", quietly = TRUE)) {
      stop("For method '", spectrum.control$method, "' you need the 'sapa' package.\n",
           "\t Please install it or user another method.")
    }    
  } else if (spectrum.control$method == "mvspec") {
    if (!requireNamespace("astsa", quietly = TRUE)) {
      stop("For method '", spectrum.control$method, "' you need the 'astsa' package.\n",
           "Please install it or user another method.")
    }
  }
  
  if (is.null(spectrum.control$smoothing)) {
    spectrum.control$smoothing <- FALSE
  }
  stopifnot(is.logical(spectrum.control$smoothing))
  
  #if (is.null(spectrum.control$taper)) {
  #  spectrum.control$taper <- 0.05
  #}
  #stopifnot(is.numeric(spectrum.control$taper),
  #          length(spectrum.control$taper) == 1,
  #          spectrum.control$taper >= 0)

  return(spectrum.control)
}