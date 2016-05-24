#' Fit Nonlinear Growth Models to Data Frame
#'
#' Determine maximum growth rates by nonlinear fits for
#' a series of experiments.
#'
#' @param formula model formula specifying dependent, independent and grouping
#'   variables in the form:
#'   \code{dependent ~ independent | group1 + group2 + \dots}.
#' @param data data frame of observational data.
#' @param time character vector with name of independent variable.
#' @param y character vector with name of dependent variable.
#' @param grouping vector of grouping variables defining subsets in the data frame.
#' @param FUN function of growth model to be fitted.
#' @param p named vector of start parameters and initial values of the growth model.
#' @param lower lower bound of the parameter vector.
#' @param upper upper bound of the parameter vector.
#' @param which vector of parameter names that are to be fitted.
#' @param method character vector specifying the optimization algorithm.
#' @param transform fit model to non-transformed or log-transformed data.
#' @param subset a specification of the rows to be used: defaults to all rows.
#' @param \dots additional parameters passed to the optimizer.
#' @param ncores number of CPU cores used for parallel computation. The number
#'   of real cores is detected automatically by default,
#'   but fort debugging purposes it could be wise to set \code{ncores = 1}.
#'   Usage of logical (hyperthreading) cores does not speed up computation.
#'
#' @return object containing the parameters of all fits.
#'
#' @family fitting functions
#'
#'
#' @examples
#'
#' data(bactgrowth)
#' splitted.data <- multisplit(value ~ time | strain + conc + replicate,
#'                  data = bactgrowth)
#'
#' ## show which experiments are in splitted.data
#' names(splitted.data)
#'
#' ## get table from single experiment
#' dat <- splitted.data[["D:0:1"]]
#'
#' fit0 <- fit_spline(dat$time, dat$value)
#'
#' fit1 <- all_splines(value ~ time | strain + conc + replicate,
#'                  data = bactgrowth, spar = 0.5)
#'
#' \donttest{
#' ## these examples require some CPU power and may take a bit longer
#'
#' ## initial parameters
#' p <- c(coef(fit0), K = max(dat$value))
#'
#' ## avoid negative parameters
#' lower = c(y0 = 0, mumax = 0, K = 0)
#'
#' ## fit all models
#' fit2 <- all_growthmodels(value ~ time | strain + conc + replicate,
#'           data = bactgrowth, FUN=grow_logistic,
#'           p = p, lower = lower, ncores = 2)
#'
#' results1 <- results(fit1)
#' results2 <- results(fit2)
#' plot(results1$mumax, results2$mumax, xlab="smooth splines", ylab="logistic")
#'
#' ## experimental: nonlinear model as part of the formula
#'
#' fit3 <- all_growthmodels(
#'           value ~ grow_logistic(time, parms) | strain + conc + replicate,
#'           data = bactgrowth, p = p, lower = lower, ncores = 2)
#'
#' ## this allows also to fit to the 'global' data set or any subsets
#' fit4 <- all_growthmodels(
#'           value ~ grow_logistic(time, parms),
#'           data = bactgrowth, p = p, lower = lower, ncores = 1)
#' plot(fit4)
#'
#' fit5 <- all_growthmodels(
#'           value ~ grow_logistic(time, parms) | strain + conc,
#'           data = bactgrowth, p = p, lower = lower, ncores = 2)
#' plot(fit5)
#' }
#'
#' @rdname all_growthmodels
#' @export
#'
all_growthmodels <- function(...) UseMethod("all_growthmodels")

#' @rdname all_growthmodels
#' @export
#'
all_growthmodels.formula <- function(formula, data,
                                     p, lower = -Inf, upper = Inf,
                                     which = names(p),
                                     FUN = NULL,
                                     method = "Marq",
                                     transform = c("none", "log"), ...,
                                     subset = NULL,
                                     ncores = detectCores(logical = FALSE)
                                     ) {

  ## FUN is now part of the formula, y ~ f(x, parms) | groups

  if (length(grep("^.*[(].*[)]", as.character(formula)[[3]]))) {   # RHS
    if (!is.null(FUN))
      warning("Nonlinear model in formula overrules value of argument FUN")

    parsed_fun <- parse_formula_nonlin(formula)

    ## simplify formula by removing a nonlinear function
    if (is.null(parsed_fun$groups)) {
      formula <- with(parsed_fun, as.formula(paste(valuevar, "~", timevar)))
    } else {
      formula <- with(parsed_fun, as.formula(paste(valuevar, "~", timevar, "|",
                                             paste(groups, collapse="+"))))
    }
    ## FUN1 = full expression; FUN2 = function name only
    FUN <- eval(parse(text = parsed_fun$FUN2))
  }

  dataset_name <- deparse(substitute(data))  # name of the dataset in the call
  X <- get_all_vars(formula, data)
  attr(X, "dataset_name") <- dataset_name

  if (!is.null(subset)) X <- X[subset, ]

  ## pass all arguments except subset
  ## grouping is formula from which y and time vars will be taken
  all_growthmodels.function(FUN=FUN, p=p, data = X, grouping = formula,
                            lower = lower, upper = upper,
                            which = which, method = method,
                            transform = transform, ...,
                            ncores = ncores)
}

#' @rdname all_growthmodels
#' @export
#'
all_growthmodels.function <-
  function(FUN, p, data, grouping = NULL, time = "time", y = "value",
                             lower = -Inf, upper = Inf,
                             which = names(p),
                             method = "Marq",
                             transform = c("none", "log"), ...,
                             ncores = detectCores(logical = FALSE)) {

  ## check arguments -----------------------------------------------------------

  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!is.character(grouping) & !inherits(grouping, "formula"))
    stop("grouping must be a formula or character vector")
  # todo:
  #if (!all(grouping %in% names(data))) stop("all grouping criteria must be column names of data")
  if (!is.function(FUN)) stop("FUN needs to be a valid growth model")
  if (!is.numericOrNull(lower) & is.numericOrNull(upper))
    stop("lower and opper must be numeric vectors or empty; lists are not possible yet")

  ## remember name of data set
  if (is.null(attr(data, "dataset_name"))) {   # inherited from former method ?
    dataset_name <- deparse(substitute(data))  # get new one
  } else {
    dataset_name <- attr(data, "dataset_name") # take old one
  }

  ## todo: consider to attach parsed formula as attr to splitted.data
  if (inherits(grouping, "formula")) {
    parsed   <- parse_formula(grouping)
    time     <- parsed$timevar
    y        <- parsed$valuevar
    grouping <- parsed$groups
  }

  ## missing groups => complete data handled as one group
  if (is.null(grouping)) {
    splitted.data <- list(data)
    names(splitted.data) <- dataset_name
    ndata <- 1
  } else {
    splitted.data <- multisplit(data, grouping)
    ndata <- length(splitted.data)
  }


  ## convert p to data frame if matrix
  if (is.matrix(p)) {
    p <- as.data.frame(p)
    ## fix empty "which" if p was a matrix
    if (is.null(which)) which <- names(p)
  }

  ## convert p to a list of the rows
  if (is.data.frame(p)) {
    p <- apply(p, 1, list)
    p <- lapply(p, unlist)
  }

  npar  <- if (is.numeric(p)) 1 else (length(p))
  if (!(npar) %in% c(1, ndata))
    stop("length of start parameters does not match number of samples")

  p1 <- if (npar == 1) p else p[[1]]
  if (!all(which %in% names(p1))) stop("parameter names from 'which' not found in p")

  nc.exist <- detectCores()
  if (ncores > nc.exist)
    warning(ncores, " cores requested but computer has only ", nc.exist)

  ## ... more checking, when necessary

  ## start of computation ------------------------------------------------------

  if (ncores == 1) {
    ## single core, p is vector with n parameter sets
    if (is.list(p)){
      fits <- lapply(seq_along(splitted.data),
             FUN = function(i) {
               fit_growthmodel(FUN, p=p[[i]],
                               time = splitted.data[[i]][,time],
                               y = splitted.data[[i]][,y],
                               lower = lower, upper = upper, #which = which,
                               method = method, transform = transform, ...
               )}
      )

    } else {
      ## single core, p is vector with 1 parameter set
      fits <- lapply(splitted.data,
        function(tmp) fit_growthmodel(FUN, p, time = tmp[,time], y = tmp[,y],
          lower = lower, upper = upper, which = which,
          method = method, transform = transform, ...
      ))
    }


  } else {
    ## multi core, p is vector with 1 or n parameter sets
    cl <- makeCluster(getOption("cl.cores", ncores))
    on.exit(stopCluster(cl))

    ## function that is to be run on the cores
    parfun <- function(X, FUNx, splitted.data = splitted.data,
                       p, time, y, lower, upper, which, method, transform, ...) {

      time <- splitted.data[[X]][ ,time]
      y    <- splitted.data[[X]][ ,y]
      p1    <- if (is.numeric(p)) p else p[[X]] # 1 or n parameter sets

      fit_growthmodel(FUNx,
                      p = p1, time = time, y = y,
                      lower = lower, upper = upper, #which = which,
                      method = method, transform = transform, ...)
    }

    ## vector of indizes of the data list
    X <- seq_along(splitted.data)

    ## multicore controller
    fits <- parLapply(cl = cl, X = X, fun = parfun,
                      FUNx = FUN,
                      splitted.data = splitted.data,
                      p=p, time = time, y = y,
                      lower = lower, upper = upper, which = which,
                      method = method, transform = transform, ...
                     )
  }
  ## names got lost during computation, so re-assign names to fits
  names(fits) <- names(splitted.data)

  ## one fit without grouping
  if (is.null(grouping)) grouping <- dataset_name

  ## create S4 object
  new("multiple_nonlinear_fits", fits = fits, grouping = grouping)
}

