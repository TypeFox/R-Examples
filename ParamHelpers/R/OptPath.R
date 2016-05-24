#' @title Create optimization path.
#'
#' @description
#' Optimizers can iteratively log their evaluated points
#' into this object. Can be converted into a data.frame with
#' \code{as.data.frame(x, discretes.as.factor = TRUE / FALSE)}.
#'
#' A optimization path has a number of path elements, where each element consists of: the value of the
#' decision variables at this point, the values of the performance measures at this point,
#' the date-of-birth (dob) of this point, the end-of-life (eol) of this point and possibly
#' an error message. See also \code{\link{addOptPathEl}}.
#'
#' For discrete parameters always the name of the value is stored as a character.
#' When you retrieve an element with \code{\link{getOptPathEl}}, this name is converted to
#' the actual discrete value.
#'
#' If parameters have associated transformation you are free to decide whether you want to
#' add x values before or after transformation, see argument \code{add.transformed.x} and
#' \code{\link{trafoOptPath}}.
#'
#' The S3 class is a list which stores at least these elements:
#' \describe{
#' \item{par.set [\code{\link{ParamSet}}]}{See argument of same name.}
#' \item{y.names [\code{character}]}{See argument of same name.}
#' \item{minimize [\code{logical}]}{See argument of same name.}
#' \item{add.transformed.x [\code{logical(1)}]}{See argument of same name.}
#' \item{env [\code{environment}]}{Environment which stores the optimization path.
#'   Contents depend on implementation.}
#' }
#'
#' @template arg_parset
#' @param y.names [\code{character}]\cr
#'   Names of performance measures that are optimized or logged.
#' @param minimize [\code{logical}]\cr
#'   Which of the performance measures in y.names should be minimized?
#'   Vector of booleans in the same order as \code{y.names}.
#' @param add.transformed.x [\code{logical(1)}]\cr
#'   If some parameters have associated transformations, are you going to
#'   add x values after they have been transformed?
#'   Default is \code{FALSE}.
#' @param include.error.message [\code{logical(1)}]\cr
#'   Should it be possible to include an error message string (or NA if no error occurred)
#'   into the path for each evaluation?
#'   This is useful if you have complex, long running objective evaluations that might fail.
#'   Default is \code{FALSE}.
#' @param include.exec.time [\code{logical(1)}]\cr
#'   Should it be possible to include execution time of evaluations
#'   into the path for each evaluation?
#'   Note that execution time could also be entered in \code{y.names} as a direct
#'   performance measure. If you use this option here, time is regarded as an extra measurement
#'   you might be curious about.
#'   Default is \code{FALSE}.
#' @param include.extra [\code{logical(1)}]\cr
#'   Should it be possible to include extra info
#'   into the path for each evaluation?
#'   Default is \code{FALSE}.
#' @name OptPath
#' @rdname OptPath
#' @family optpath
NULL

makeOptPath = function(par.set, y.names, minimize, add.transformed.x = FALSE,
  include.error.message = FALSE, include.exec.time = FALSE, include.extra = FALSE) {

  n.y = length(y.names)
  ok = c("numeric", "integer", "numericvector", "integervector", "logical",
    "logicalvector", "discrete", "discretevector", "character", "charactervector")
  if(length(par.set$pars) > length(filterParams(par.set, type = ok)$pars))
    stop("OptPath can currently only be used for: ", paste(ok, collapse = ","))
  x.names = getParamIds(par.set)
  # be really sure that x and y columns are uniquely named
  x.names2 = c(getParamIds(par.set, with.nr = TRUE), getParamIds(par.set, with.nr = FALSE))
  if (length(intersect(x.names2, y.names)) > 0)
    stop("'x.names' and 'y.names' must not contain common elements!")
  if (length(minimize) != n.y)
    stop("'y.names' and 'minimize' must be of the same length!")
  if (is.character(names(minimize)) && !setequal(names(minimize), y.names))
    stop("Given names for 'minimize' must be the same as 'y.names'!")
  if (is.null(names(minimize)))
    names(minimize) = y.names
  if (any(c("dob", "eol", "error.message") %in% (union(x.names, y.names))))
    stop("'dob', 'eol' and 'error.message' are not allowed in parameter names or 'y.names'!")
  ee = new.env()
  ee$dob = ee$eol = integer(0)

  # potentially init error.message and exec.time in env
  ee$error.message = if (include.error.message) character(0L) else NULL
  ee$exec.time = if (include.exec.time) numeric(0L) else NULL
  ee$extra = if (include.extra) list() else NULL

  makeS3Obj("OptPath",
    par.set = par.set,
    y.names = y.names,
    minimize = minimize,
    add.transformed.x = add.transformed.x,
    env = ee
  )
}

#' @export
print.OptPath = function(x, ...) {
  n = getOptPathLength(x)
  em = x$env$error.message
  et = x$env$exec.time
  ex = x$env$extra
  catf("Optimization path")
  catf("  Dimensions: x = %i/%i, y = %i",
    length(x$par.set$pars), sum(getParamLengths(x$par.set)), length(x$y.names))
  catf("  Length: %i", n)
  catf("  Add x values transformed: %s", x$add.transformed.x)
  s = if (is.null(em)) ""  else sprintf(" Errors: %i / %i.", sum(!is.na(em)), n)
  catf("  Error messages: %s.%s", !is.null(em), s)
  s = if (is.null(et)) {
    s = ""
  } else {
    ntimes = sum(!is.na(et))
    ntime.nas = sum(is.na(et))
    # no non-na exec times in path
    if (ntimes == 0L) {
      et1 = 0
      et2 = 0
     } else {
      et1 = min(et, na.rm = TRUE)
      et2 = max(et, na.rm = TRUE)
    }
    s = sprintf(" Range: %g - %g. %i NAs.", et1, et2, ntime.nas)
  }
  catf("  Exec times: %s.%s", !is.null(et), s)
  if (!is.null(ex))
  catf("  Extras: %i columns", ifelse(length(ex) > 0L, length(ex[[1]]), NA))
}



