#' @name checking
#' 
#' @param x            scalar or vector of quantiles
#' @param prob         scalar or vector of probability  
#' @param n            scalar sample size
#' @param param        scalar or vector of parameters
#' @param phiu         scalar or vector of phiu (logical, NULL or 0-1 exclusive)
#' @param ns           vector of lengths of parameter vectors
#' @param nparam       acceptable length of (non-scalar) vectors of parameter vectors
#' @param logicarg     logical input argument
#' @param textarg      character input argument
#' @param inputn       vector of input lengths
#' @param allowvec     logical, where TRUE permits vector
#' @param allownull    logical, where TRUE permits NULL values
#' @param allowmiss    logical, where TRUE permits missing input
#' @param allowna      logical, where TRUE permits NA and NaN values
#' @param allowzero    logical, where TRUE permits zero values (positive vs non-negative)
#' @param allowscalar  logical, where TRUE permits scalar (as opposed to vector) values
#' @param allowinf     logical, where TRUE permits +/-Inf values
#' @param allowfalse   logical, where TRUE permits FALSE (and TRUE) values
#' @param method       optimisation method (see \code{\link[stats:optim]{optim}})
#' @param control      optimisation control list (see \code{\link[stats:optim]{optim}})
#' @param bcmethod     boundary correction method
#' @param nn           non-negativity correction method (simple boundary correction only)
#' @param offset       offset added to kernel centres (logtrans only) or \code{NULL}
#' @param beta         vector of B-spline coefficients (required)
#' @param xrange       vector of minimum and maximum of B-spline (support of density)
#' @param nseg         number of segments between knots
#' @param degree       degree of B-splines (0 is constant, 1 is linear, etc.)
#' @param design.knots spline knots for splineDesign function
#' 
#' @title Internal functions for checking function input arguments
#'
#' @description Functions for checking the input arguments to functions, so that main functions
#' are more concise. They will stop when an inappropriate input is found.
#' 
#' These function are visible and operable by the user. But they should be used with caution, as no
#' checks on the input validity are carried out.
#' 
#' For likelihood functions you will often not want to stop on finding a non-positive values for
#' positive parameters, in such cases use \code{\link[evmix:checking]{check.param}} rather than 
#' \code{\link[evmix:checking]{check.posparam}}.
#' 
#' @return The checking functions will stop on errors and return no value. The only exception is
#' the \code{\link[evmix:checking]{check.inputn}} which outputs the maximum vector length.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}.
#'
#' @aliases checking check.quant check.prob check.n check.param check.posparam
#'  check.nparam check.logic check.text check.inputn
#' @family checking
#' 
NULL

#' @export
#' @rdname checking
check.param <- function(param, allowvec = FALSE, allownull = FALSE, allowmiss = FALSE,
                        allowna = FALSE, allowinf = FALSE) {
  
  fcall = match.call()
  
  if (missing(param) & !allowmiss) stop(paste("missing parameter input from call", deparse(fcall)))
  
  # If non-missing then need to check parameter value
  if (!missing(param)) {
    
    param.name = deparse(fcall[[2]])
    fcall = paste(deparse(fcall), collapse = "")
    
    if (is.null(param) & !allownull)
      stop(paste("parameter input", ifelse(is.null(param.name), "(unnamed)", param.name),
                 "is NULL from call", fcall))

    if (!is.null(param)) {
            
      if (length(param) == 0)
        stop(paste("parameter input", param.name, "is empty from call", fcall))
    
      # Special case - if all elements are NA the vector/scalar can be logical
      if (all(is.na(param))) {
        if (!allowna) stop(paste("NA/NaN values for parameter input", param.name, "not permitted from call", fcall))
      } else {
        if (mode(param) != "numeric")
          stop(paste("parameter input", param.name, "must be a non-empty numeric",
                     ifelse(allowvec, "vector", "scalar"), "from call", fcall))
      
        if (any(is.na(param)) & !allowna)
          stop(paste("NA/NaN values for parameter input", param.name, "not permitted from call", fcall))
  
        if (any(is.infinite(param)) & !allowinf)
          stop(paste("infinite values for parameter input", param.name, "not permitted from call", fcall))
      }
      if ((length(param) > 1) & !allowvec)
        stop(paste("parameter input", param.name, "must be scalar from call", fcall))
    }
  }
}

#' @export
#' @rdname checking
check.posparam <- function(param, allowvec = FALSE, allownull = FALSE, allowmiss = FALSE,
                        allowna = FALSE, allowinf = FALSE, allowzero = FALSE) {

  check.param(param, allowvec, allownull, allowmiss, allowna, allowinf)
    
  if (!missing(param)) {
    
    if (!is.null(param)) {
      if (!all(is.na(param))) {
        fcall = match.call()
        param.name = deparse(fcall[[2]])
        fcall = paste(deparse(fcall), collapse = "")
      
        if (any(param < 0, na.rm = TRUE)) 
          stop(paste("non-positive parameter", param.name, "not permitted from call", fcall))
      
        if (any(param == 0, na.rm = TRUE) & !allowzero) 
          stop(paste("zero parameter", param.name, "not permitted from call", fcall))
      }
    }
  }
}

#' @export
#' @rdname checking
check.quant <- function(x, allownull = FALSE, allowna = FALSE, allowinf = FALSE) {
  # quantiles vector, must be non-missing and non-empty, does not check if vector/scalar
  
  fcall = match.call()
  
  if (missing(x)) stop(paste("missing quantiles from call", deparse(fcall)))
  
  x.name = deparse(fcall[[2]])
  fcall = paste(deparse(fcall), collapse = "")
  
  if (is.null(x) & !allownull)
    stop(paste("quantiles input", ifelse(is.null(x.name), "(unnamed)", x.name), "is NULL from call", fcall))
  
  if (!is.null(x)) {
    
    if (length(x) == 0)
      stop(paste("quantiles input", x.name, "is empty from call", fcall))

      # Special case - if all elements are NA the vector/scalar can be logical
    if (all(is.na(x))) {
      if (!allowna) stop(paste("NA/NaN values for quantiles input", x.name, "not permitted from call", fcall))
    } else {
      if (mode(x) != "numeric")
        stop(paste("quantiles input", x.name, "must be a non-empty numeric vector from call", fcall))
      
      if (any(is.na(x)) & !allowna)
        stop(paste("NA/NaN values for quantiles input", x.name, "not permitted from call", fcall))
      
      if (any(is.infinite(x)) & !allowinf)
        stop(paste("infinite values for quantiles input", x.name, "not permitted from call", fcall))
    }
  }
}

#' @export
#' @rdname checking
check.prob <- function(prob, allownull = FALSE, allowna = FALSE) {
  # probability vector, must be non-missing and non-empty, does not check if vector/scalar

  check.quant(prob, allownull, allowna)
  
  fcall = match.call()
  prob.name = deparse(fcall[[2]])
  fcall = paste(deparse(fcall), collapse = "")
  
  if (!is.null(prob)) {
    if (any(prob < 0, na.rm = TRUE) | any(prob > 1, na.rm = TRUE))
      stop(paste("probability input", prob.name, "must be between 0 and 1 (inclusive) from call", fcall))
  }
}

#' @export
#' @rdname checking
check.n <- function(n, allowzero = FALSE) {
  # positive integer only permitted by default
  # non-negative when allowzero = TRUE

  fcall = match.call()

  if (missing(n)) stop("sample size input is missing")
  
  n.name = deparse(fcall[[2]])
  fcall = paste(deparse(fcall), collapse = "")
    
  if (is.null(n))
    stop(paste("sample size input", ifelse(is.null(n.name), "(unnamed)", n.name),
               "is NULL from call", fcall))

  if (length(n) == 0)
    stop(paste("sample size input", n.name, "is empty from call", fcall))
    
  if (length(n) != 1)
    stop(paste("sample size input", n.name, "must be scalar",
               ifelse(allowzero, "non-negative", "positive"), "integer from call", fcall))
    
  if (mode(n) != "numeric")
    stop(paste("sample size input", n.name, "must be a scalar",
               ifelse(allowzero, "non-negative", "positive"), "integer from call", fcall))
  
  if (!is.finite(n))
    stop(paste("sample size input", n.name, "is not finite from call", fcall))
  
  if (abs(n - round(n)) > sqrt(.Machine$double.eps))
    stop(paste("sample size input", n.name, "must be a scalar",
               ifelse(allowzero, "non-negative", "positive"), "integer from call", fcall))

  if (n < 0)
    stop(paste("sample size input", n.name, "must be a scalar",
               ifelse(allowzero, "non-negative", "positive"), "integer from call", fcall))
  
  if ((n < 1) & !allowzero)
    stop(paste("sample size input", n.name, "must be a scalar positive integer from call", fcall))
}

#' @export
#' @rdname checking
check.logic <- function(logicarg, allowvec = FALSE, allowna = FALSE) {
  # only logical, numeric(0)/NULL/NaN not permitted as non-logical
  
  fcall = match.call()
  
  if (missing(logicarg)) stop(paste("missing logical input from call", deparse(fcall)))

  logic.name = deparse(fcall[[2]])
  fcall = paste(deparse(fcall), collapse = "")
  
  if (is.null(logicarg)) # NULL not permitted
    stop(paste("NULL value for input", logic.name, "not permitted, must be logical",
               ifelse(allowvec, "vector", "scalar"), "from call", fcall))

  if (length(logicarg) == 0) # numeric(0)
    stop(paste("empty value for input", logic.name, "not permitted, must be logical", 
               ifelse(allowvec, "vector", "scalar"), "from call", fcall))
    
  if (any(is.nan(logicarg))) # NaN not permitted
    stop(paste("NaN value for input", logic.name, "not permitted, must be logical",
               ifelse(allowvec, "vector", "scalar"), "from call", fcall))
  
  if (any(is.infinite(logicarg))) # infinite values not permitted
    stop(paste("infinite value for input", logic.name, "not permitted, must be logical",
               ifelse(allowvec, "vector", "scalar"), "from call", fcall))

  # above 3 checks would also be caught by following, but above is more helpful to user
  if (mode(logicarg) != "logical") # non-logical not permitted
    stop(paste("non-logical input", logic.name, "not permitted, must be logical",
               ifelse(allowvec, "vector", "scalar"), "from call", fcall))

  if (any(is.na(logicarg)) & !allowna)
    stop(paste("NA values for logical input", logic.name, "not permitted from call", fcall))

  if ((length(logicarg) > 1) & !allowvec)
    stop(paste("input", logic.name, "must be logical scalar from call", fcall))
}

#' @export
#' @rdname checking
check.nparam <- function(ns, nparam = 1, allownull = FALSE, allowmiss = FALSE) {

  fcall = match.call()
  
  if (missing(ns) & !allowmiss) stop(paste("missing parameter input from call", deparse(fcall)))
  
  # If non-missing then need to check parameter value
  if (!missing(ns)) {
    
    ns.name = deparse(fcall[[2]])
    fcall = paste(deparse(fcall), collapse = "")

    check.n(nparam, allowzero = TRUE)
    
    if (is.null(ns) & !allownull)
      stop(paste("parameter input", ifelse(is.null(ns.name), "(unnamed)", ns.name),
                 "is NULL from call", fcall))
    
    if (length(ns) == 0)  { # empty/NULL OK if nparam == 0
      if (!is.null(ns) & (nparam != 0))
        stop(paste("parameter input", ns.name, "is empty from call", fcall))      
    } else {

      # Special case - if all elements are NA the vector/scalar can be logical
      if (!all(is.na(ns))) {
        if (mode(ns) != "numeric")
          stop(paste("parameter input", ns.name, "must be numeric",
                     ifelse(nparam > 1, "vector", "scalar"), "from call", fcall))
      }
      if (length(ns) != nparam)
        stop(paste("parameter input", ns.name, "must be numeric",
                     ifelse(nparam > 1, "vector", "scalar"), "of length", nparam, "from call", fcall))
    }
  }
}

#' @export
#' @rdname checking
check.inputn <- function(inputn, allowscalar = FALSE, allowzero = FALSE) {
  # parameter lengths consistency check, mixture of scalar/vectors permitted by default
  # vector of parameter lengths (scalar irrelevant) in main input

  fcall = match.call()
  
  if (missing(inputn)) stop(paste("missing vector lengths (inputn) from call", deparse(fcall)))
      
  inputn.name = deparse(fcall[[2]])
  fcall = paste(deparse(fcall), collapse = "")

  if (is.null(inputn))
    stop(paste("vector lengths", ifelse(is.null(inputn.name), "(unnamed)", inputn.name),
               "is NULL from call", fcall))
  
  if (length(inputn) == 0)
    stop(paste("vector lengths", inputn.name, "is empty vector from call", fcall))
          
  if (mode(inputn) != "numeric")
    stop(paste("vector lengths", inputn.name, "must be",
      ifelse(allowzero, "non-negative", "positive"), "integer vector from call", fcall))
          
  if (any(!is.finite(inputn)))
    stop(paste("vector lengths", inputn.name, "must be",
      ifelse(allowzero, "non-negative", "positive"), "integer vector from call", fcall))

  if (any(abs(inputn - round(inputn)) > sqrt(.Machine$double.eps)))
    stop(paste("vector lengths", inputn.name, "must be",
      ifelse(allowzero, "non-negative", "positive"), "integer vector from call", fcall))
  
  if (any(inputn < 0))
    stop(paste("vector lengths", inputn.name, "must be",
      ifelse(allowzero, "non-negative", "positive"), "integer vector from call", fcall))

  if (any(inputn == 0) & (!allowzero))
    stop(paste("vector lengths", inputn.name, "must be positive integer vector from call", fcall))

  n = max(inputn)

  if (!isTRUE(all.equal(inputn, rep(n, length(inputn))))) {
    if (isTRUE(all.equal(inputn[inputn != 1], rep(n, length(inputn[inputn != 1]))))) {
      if (!allowscalar) 
        stop(paste("vector lengths", inputn.name, "must be",
                   ifelse(allowzero, "non-negative", "positive"),
                   "integer vector (with all vectors of same length) from call", fcall))
    } else {
      stop(paste("vector lengths", inputn.name, "must be",
                 ifelse(allowzero, "non-negative", "positive"),
                 "integer vector (with all vectors of same length, not scalars) from call", fcall))  
    }
  }

  return(n)
}

#' @export
#' @rdname checking
check.text <- function(textarg, allowvec = FALSE, allownull = FALSE) {

  fcall = match.call()
  
  if (missing(textarg)) stop(paste("missing character input (textarg) from call", deparse(fcall)))

  param.name = deparse(fcall[[2]])
  fcall = paste(deparse(fcall), collapse = "")

  if (is.null(textarg) & !allownull)
    stop(paste("character input", ifelse(is.null(param.name), "(unnamed)", param.name),
               "is NULL from call", fcall))

  if (!is.null(textarg)) {
    
    if (length(textarg) == 0)
      stop(paste("character input", ifelse(is.null(param.name), "(unnamed)", param.name),
                   "is empty from call", fcall))

    if (any(is.na(textarg))) {
      stop(paste("NA/NaN values for character input", param.name, "not permitted from call", fcall))
    } else {
      if (mode(textarg) != "character")
        stop(paste("non-character values for character input", param.name, "not permitted from call", fcall))
    
      if ((length(textarg) > 1) & !allowvec)
        stop(paste("character input", param.name, "must be single string from call", fcall))
    }
  }
}

#' @export
#' @rdname checking
check.phiu <- function(phiu, allowvec = FALSE, allownull = FALSE, allowfalse = FALSE) {
  # vector or scalar of phiu, must be non-missing and non-empty
  # NA/NaN and infinite values never permitted
  # TRUE/FALSE or 0-1 permitted in likelihood, as full sample data is available
  # TRUE or 0-1 (exclusive) in other functions, as sample proportion (FALSE) not available

  fcall = match.call()
  
  if (missing(phiu)) stop(paste("missing phiu input from call", deparse(fcall)))
  
  phiu.name = deparse(fcall[[2]])
  fcall = paste(deparse(fcall), collapse = "")
  
  if (allowfalse) {
    phiuerror = "phiu must be either TRUE for bulk model based tail fraction approach, 
      or between 0 and 1 (exclusive) when using parameterised tail fraction approach"
  } else {
    phiuerror = "phiu must be either TRUE for bulk model based tail fraction approach, 
      FALSE for parameterised tail fraction approach (estimated using sample proportion)
      or between 0 and 1 (exclusive) for fixed tail fraction"
  }

  if (is.null(phiu) & !allownull)
    stop(paste("phiu input", ifelse(is.null(phiu.name), "(unnamed)", phiu.name),
               "is NULL from call", fcall, ".", phiuerror))
  
  if (!is.null(phiu)) {
    
    if (length(phiu) == 0)
      stop(paste("phiu input", phiu.name, "is empty from call", fcall, ".", phiuerror))
  
    if ((mode(phiu) != "numeric") & (mode(phiu) != "logical"))
        stop(paste("phiu input", phiu.name, "must be a non-empty numeric or logical",
                   ifelse(allowvec, "vector", "scalar"), "from call", fcall, ".", phiuerror))
        
    if (any(is.na(phiu)))
      stop(paste("NA/NaN values for phiu input", phiu.name, "not permitted from call", fcall, ".", phiuerror))

    if (any(is.infinite(phiu)))
      stop(paste("infinite values for phiu input", phiu.name, "not permitted from call", fcall, ".", phiuerror))
  
    if ((length(phiu) > 1) & !allowvec)
      stop(paste("phiu input", phiu.name, "must be scalar from call", fcall, ".", phiuerror))

    if (is.logical(phiu)) {
      if(!allowfalse & any(!phiu))
        stop(paste("FALSE values for phiu input", phiu.name, "not permitted from call", fcall, ".", phiuerror))
    } else {
      check.prob(phiu)
    }
  }
}

#' @export
#' @rdname checking
check.optim <- function(method) {
  # checks optim method specification

  check.text(method)
  
  allmethod = eval(formals(optim)$method)
  if (!(method %in% allmethod))
    stop(paste(c("optim method must be one of", allmethod), collapse = " "))
}

#' @export
#' @rdname checking
check.control <- function(control) {
  # checks optim method specification

  if (missing(control)) stop("optim control list is missing")

  if (mode(control) != "list")
    stop("optim control must be list, see optim function for details")

  allcontrol = c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol",
    "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "tmax", "temp")  
  if (any(!(names(control) %in% allcontrol)))
    stop(paste(c("optim control list must be containly only", allcontrol), collapse = " "))
}

#' @export
#' @rdname checking
check.bcmethod <- function(bcmethod) {
  # checks boundary correction method specification

  check.text(bcmethod)

  allbcmethod = c("simple", "cutnorm", "renorm", "reflect", "logtrans", 
    "beta1", "beta2", "gamma1", "gamma2", "copula")
  if (!(bcmethod %in% allbcmethod))
    stop(paste(c("boundary correction method must be one of", allbcmethod), collapse = " "))
}

#' @export
#' @rdname checking
check.nn <- function(nn) {
  # checks non-negative correction method specification

  check.text(nn)

  allnn = c("none", "zero", "jf96")
  if (!(nn %in% allnn))
    stop(paste(c("non-negative correction method must be one of", allnn), collapse = " "))
}

#' @export
#' @rdname checking
check.offset <- function(offset, bcmethod, allowzero = FALSE) {
  # checks offset for different bcmethods

  check.posparam(offset, allownull = TRUE)
  check.bcmethod(bcmethod)

  if (bcmethod != "logtrans") {
    if (!is.null(offset)) warning("offset only relevant for logtrans method")
  } else {
    if (is.null(offset)) stop("offset must be provided for logtrans method")
    
    if (offset < 0) stop(paste("offset must be", ifelse(allowzero, "non-negative", "positive")))
    
    if ((offset == 0) & !allowzero ) stop("offset must be positive")
  }
}

#' @export
#' @rdname checking
check.design.knots <- function(beta, xrange, nseg, degree, design.knots) {
  # checks knot specification and calculate them if not pre-specified
  
  # Outputs a list containing the design.knots and corrected inputs
  
  # Two options for specifying knots:
  #    1) design.knots vector
  #    2) xrange, nseg and degree
  # if both provided then design.knots is used
  if (is.null(design.knots) & is.null(xrange)) stop("Either design.knots or xrange must be specified")

  if (is.null(design.knots)) {
    # check x-range
    if (length(xrange) != 2) stop("knot range in xrange must be vector of length 2")
    if (diff(xrange) <= 0) stop("knot range in xrange must have positive width")
    
    # consistent with Eilers and Marx the "P-spline masters":
    # defaults to regular knots and not natural B-splines,
    # so each B-spline is just shifted/spliced version of each other
    dx = diff(xrange)/nseg # regular knot spacing
    design.knots = seq(xrange[1] - degree * dx, xrange[2] + degree * dx, by = dx)
    
  } else {
    # if knots specified, they must be sorted
    if (is.unsorted(design.knots)) {
      design.knots = sort(design.knots)
    } else {
      if (design.knots[1] > design.knots[length(design.knots)])
        design.knots = rev(design.knots)
    }
    
    # degree determinable by difference between design.knots and coefficients
    if ((length(design.knots) - length(beta) - 1) != degree) {
      warning(paste("Degree ", degree, 
                    " is inconsistent with difference in coefficients and design knots #(design.knots) - #(beta) - 1 = ",
                    length(design.knots) - length(beta) - 1, ", so degree is reset", sep=""))
      degree = length(design.knots) - length(beta) - 1
    }

    # number of segments degree determinable from design.knots and degree
    if (nseg != (length(design.knots) - 1 - degree*2)) {
      warning(paste("Number of segments ", nseg, " is inconsistent with design knots where length(design.knots) - 1 - degree*2= ",
                    length(design.knots) - 1 - degree*2, ", so nseg reset", sep=""))
      nseg = length(design.knots) - 1 - degree*2
    }

    # xrange also determined by design knots and degree
    if (!is.null(xrange)) {
      if (!isTRUE(all.equal(xrange, design.knots[c(degree + 1, length(design.knots) - degree)]))) {
        warning(paste("Interpolation range (", xrange[1], ", ", xrange[2], ") is inconsistent with design knots (",
                      design.knots[degree + 1], ", " , design.knots[length(design.knots) - degree],
                      "), so xrange reset", sep=""))
        xrange = design.knots[c(degree + 1, length(design.knots) - degree)]
      }
    } else {
        xrange = design.knots[c(degree + 1, length(design.knots) - degree)]
    }
  }
  
  if (length(beta) != (nseg + degree)) 
    stop(paste("Number of coefficients ", length(beta), " is not equal to number of B-splines nseg + degree = ", nseg + degree, sep=""))

  return(list(xrange = xrange, nseg = nseg, degree = degree, design.knots = design.knots))
}
