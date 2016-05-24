#' Fit Exponential Growth Model with Smoothing Spline
#'
#' Determine maximum growth rates from log-linear part of the growth curve for
#' a series of experiments by using smoothing splines.
#'
#' @param formula model formula specifying dependent, independent and grouping
#'   variables in the form:
#'   \code{dependent ~ independent | group1 + group2 + \dots}.
#' @param data data frame of observational data.
#' @param grouping vector of grouping variables defining subsets in the data frame.
#' @param time character vectors with name independent variable.
#' @param y character vector with name of dependent variable.
#' @param optgrid number of steps on the x-axis used for searching the maximum
#'  of the first derivative of the spline.
#'  The default should work in most cases, as long as the data are equally spaced.
#'  A smaller number may lead to non-detectable speed-up, but has the risk that
#'  the search is trapped in a local minimum.
#' @param subset a specification of the rows to be used: defaults to all rows.
#' @param \dots other parameters passed to \code{\link{smooth.spline}}, see details.
#'
#' @return object with parameters of the fit.
#'
#' @details The method was inspired by an algorithm of Kahm et al. (2010),
#'   with different settings and assumptions. In the moment, spline fitting
#'   is always done with log-transformed data, assuming exponential growth
#'   at the time point of the maximum of its first derivative.
#'
#'   All the hard work is done by function \code{\link{smooth.spline}} from package
#'   \pkg{stats}, that is highly user configurable. Normally, smoothness is
#'   automatically determined via cross-validation. This works well in many cases,
#'   whereas manual adjustment is required otherwise, e.g. by setting \code{spar}
#'   to a fixed value \eqn{[0,1]} that also disables cross-validation.
#'   A typical case where cross validation does not work is, if time dependent
#'   measurements are taken as pseudoreplicates from the same experimental unit.
#'
#' @references
#'
#' Kahm, M., Hasenbrink, G., Lichtenberg-Frate, H., Ludwig, J., Kschischo, M.
#' 2010. grofit: Fitting Biological Growth Curves with R.
#' Journal of Statistical Software, 33(7), 1-21. URL
#' \url{http://www.jstatsoft.org/v33/i07/}
#'
#' @family fitting functions
#'
#' @examples
#'
#' data(bactgrowth)
#' L <- all_splines(value ~ time | strain + conc + replicate,
#'                  data = bactgrowth, spar = 0.5)
#'
#' par(mfrow=c(4, 3))
#' plot(L)
#' results <- results(L)
#' xyplot(mumax ~ log(conc + 1)|strain, data=results)
#'
#' ## fit splines at lower grouping levels
#' L2 <- all_splines(value ~ time | conc + strain,
#'                     data = bactgrowth, spar = 0.5)
#' plot(L2)
#'
#' ## total data set without any grouping
#' L3 <- all_splines(value ~ time,
#'                     data = bactgrowth, spar = 0.5)
#' par(mfrow=c(1, 1))
#' plot(L3)
#'
#' @rdname all_splines
#' @export
#'
all_splines <- function(...) UseMethod("all_splines")

#' @rdname all_splines
#' @export
#'
all_splines.formula <- function(formula, data=NULL, optgrid = 50, subset=NULL,  ...) {

  dataset_name <- deparse(substitute(data))  # name of the dataset in the call
  X <- get_all_vars(formula, data)
  attr(X, "dataset_name") <- dataset_name

  if (!is.null(subset)) X <- X[subset, ]
  all_splines.data.frame(data = X, grouping = formula, optgrid = optgrid, ...)
}

#' @rdname all_splines
#' @export
#'
all_splines.data.frame <-
  function(data, grouping=NULL, time = "time", y = "value",  optgrid = 50, ...) {

    ## remember name of data set
    if (is.null(attr(data, "dataset_name"))) {   # inherited from former method ?
      dataset_name <- deparse(substitute(data))  # get new one
    } else {
      dataset_name <- attr(data, "dataset_name") # take old one
    }

    ## check and parse grouping if formula
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

    #splitted.data <- multisplit(data, grouping)

    ## supress warnings, esp. in case of "perfect fit"
    fits <- lapply(splitted.data,
                   function(tmp)
                     suppressWarnings(fit_spline(tmp[,time], tmp[,y],
                                                 optgrid = optgrid, ...)))

    ## one fit without grouping
    if (is.null(grouping)) grouping <- dataset_name

    new("multiple_smooth.spline_fits", fits = fits, grouping = grouping)
  }
