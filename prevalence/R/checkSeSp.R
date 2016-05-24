checkSeSp <-
function(x, type = "prob") {
  distr <- NULL
  param <- NULL

  if (length(x) > 6)
    warning("Too many parameters specified")

  if (length(x) == 1) {
    distr <- "fixed"

    ## probability is constrained in (0,1)
    if (type == "prob" && (x < 0 | x > 1))
      stop("Parameter 'par' of fixed distribution must be",
           " numeric value between 0 and 1")

    ## covariance is constrained in +/- 2^-h
    if (is.numeric(type) && (abs(x) > 2^-type))
      stop("Parameter 'par' of fixed distribution must be",
           " numeric value between +/-", 2^-type)

    param <- x

  } else {
    distr <- x$dist
  }

  if (is.list(x) && length(unlist(x)) > length(x))
    stop("Parameters cannot be specified as vectors")

  if (is.null(distr))
    stop("No distribution specified")
  if (!is.character(distr))
    stop("Invalid distribution specified")
  distr <- tolower(distr)
  if (!any(c("fixed", "uniform", "pert", "beta", "beta-expert") == distr))
    stop("Distribution must be",
         " 'fixed', 'uniform', 'pert', 'beta' or 'beta-expert'")

  ## Fixed distribution
  if (distr == "fixed" & is.null(param)) {
    if (length(x) > 2)
      warning("A fixed distribution requires only 1 parameter")
    if (is.null(x$par))
      stop("Parameter 'par' not specified")

    ## probability is constrained in (0,1)
    if (type == "prob" && (x$par < 0 | x$par > 1))
      stop("Parameter 'par' of fixed distribution must be",
           " numeric value between 0 and 1")

    ## covariance is constrained in +/- 2^-h
    if (is.numeric(type) && (abs(x$par) > 2^-type))
      stop("Parameter 'par' of fixed distribution must be",
           " numeric value between +/-", 2^-type)

    param <- x$par
  }

  ## Uniform distribution
  if (distr == "uniform") {
    if (length(x) > 3)
      warning("A uniform distribution requires only 2 parameters")
    if (length(x) < 3)
      warning("A uniform distribution requires 2 parameters")
    if (is.null(x$min))
      stop("Parameter 'min' not specified")
    if (is.null(x$max))
      stop("Parameter 'max' not specified")
    if (x$min > x$max)
      stop("'min' of uniform distribution cannot be larger than 'max'")

    ## probability is constrained in (0,1)
    if (type == "prob" && (any(c(x$min, x$max) < 0) | any(c(x$min, x$max) > 1)))
      stop("Parameters of uniform distribution must be",
           " numeric values between 0 and 1")

    ## covariance is constrained in +/- 2^-h
    if (is.numeric(type) && any(abs(c(x$min, x$max)) > 2^-type))
      stop("Parameters of uniform distribution must be",
           " numeric values between +/-", 2^-type)

    param <- c(x$min, x$max)
  }

  ## Beta distribution
  if (distr == "beta") {
    ## covariance is constrained in +/- 2^-h
    if (is.numeric(type))
      stop("Beta distribution not allowed for covariance parameters")

    if (length(x) > 3)
      warning("A beta distribution requires only 2 parameters")
    if (length(x) < 3)
      warning("A beta distribution requires 2 parameters")
    if (is.null(x$alpha))
      stop("Parameter 'alpha' not specified")
    if (is.null(x$beta))
      stop("Parameter 'beta' not specified")
    if (any(c(x$alpha, x$beta) <= 0))
      stop("Parameters of beta distribution must be",
           " numeric values larger than 0")

    param <- c(x$alpha, x$beta)
  }

  ## Beta-Expert distribution
  if (distr == "beta-expert") {
    ## covariance is constrained in +/- 2^-h
    if (is.numeric(type))
      stop("Beta-Expert distribution not allowed for covariance parameters")

    if (is.null(x$mode) & is.null(x$mean))
      stop("At least 'mode' or 'mean' must be specified")
    if (!is.null(x$mode) & !is.null(x$mean))
      stop("'mode' and 'mean' cannot both be specified")
    method <- c("mode", "mean")[c(!is.null(x$mode), !is.null(x$mean))]
    best <- ifelse(method == "mode", x$mode, x$mean)
    if (is.null(x$lower) & is.null(x$upper))
      stop("At least 'lower' or 'upper' must be specified")
    if (is.null(x$p))
      stop("Parameter 'p' not specified")
    target <- c(x$lower, x$upper)[c(!is.null(x$lower), !is.null(x$upper))]

    ## probability is constrained in (0,1)
    if (type == "prob" && (any(c(best, x$p, x$target) < 0) | any(c(best, x$p, x$target) > 1)))
      stop("Parameters of beta-expert distribution must be",
           " numeric values between 0 and 1")

    if (!is.null(x$lower))
      if (x$lower > x$m)
        stop("'lower' cannot be larger than 'm'")
    if (!is.null(x$upper))
      if (x$upper < x$m)
        stop("'upper' cannot be smaller than 'm'")
    if (!is.null(x$lower) & !is.null(x$upper))
      if (x$lower > x$upper)
        stop("'lower' cannot be larger than 'upper'")

    distr <- "beta"
    if (is.null(x$upper)) {
      param <-
        betaExpert(best = best, method = method, lower = x$lower, p = x$p)
    } else if (is.null(x$lower)) {
      param <-
        betaExpert(best = best, method = method, upper = x$upper, p = x$p)
    } else {
      param <-
        betaExpert(best = best, method = method,
                   lower = x$lower, upper = x$upper, p = x$p)
    }
  }

  ## Beta-PERT distribution
  if (distr == "pert") {
    if (length(x) > 6)
      warning("A PERT distribution requires maximum 5 parameters")
    if (length(x) < 4)
      warning("A PERT distribution requires at least 3 parameters")
    if (is.null(x$a))
      stop("Parameter 'a' not specified")
    if (is.null(x$m))
      stop("Parameter 'm' not specified")
    if (is.null(x$b))
      stop("Parameter 'b' not specified")

    ## probability is constrained in (0,1)
    if (type == "prob" && (any(c(x$a, x$m, x$b) < 0) | any(c(x$a, x$m, x$b) > 1)))
      stop("Parameters of PERT distribution must be",
           " numeric values between 0 and 1")

    ## covariance is constrained in +/- 2^-h
    if (is.numeric(type) && any(abs(c(x$a, x$m, x$b)) > 2^-type))
      stop("Parameters of PERT distribution must be",
           " numeric values between +/-", 2^-type)

    if (x$a > x$m)
      stop("'a' of PERT distribution cannot be larger than 'm'")
    if (x$m > x$b)
      stop("'m' of PERT distribution cannot be larger than 'b'")
    pertK <- ifelse(is.null(x$k), 4, x$k)
    pertM <- ifelse(is.null(x$method), "classic", x$method)
    param <- c(x$a, x$m, x$b, pertK)
    distr <- c(distr, pertM)
  }

  return(list(d = distr, p = as.numeric(param)))
}