#' install.dates
#'
#' @description Download and install the latest versions of packages hosted on
#'  CRAN as of a specific date from the MRAN server.
#'
#' @param pkgs character vector of the names of packages that should be
#'  downloaded and installed
#'
#' @param dates character or Date vector of the dates for which to install the
#' latest versions of \code{pkgs}. If a data vector, it must be in the format
#' 'yyyy-mm-dd', e.g. '2014-09-17'. If this has the same length as \code{pkgs}
#'  versions will correspond to those packages. If this has length one
#'  the same version will be used for all packages. If it has any other
#'  length an error will be thrown. Dates before 2014-09-17 will cause an error
#'  as MRAN does not archive before that date.
#'
#' @param lib character vector giving the library directories where to
#'  install the packages. Recycled as needed. If missing, defaults to the
#'  first element of \code{\link{.libPaths}()}.
#'
#' @param \dots other arguments to be passed to \code{\link{install.packages}}.
#'  The arguments \code{repos} and \code{contriburl} (at least) will
#'  be ignored as the function uses the MRAN server to retrieve package versions.
#'
#' @export
#' @name install.dates
#' @examples
#'
#' \dontrun{
#'
#' # install yesterday's version of checkpoint
#' install.dates('checkpoint', Sys.Date() - 1)
#'
#' # install yesterday's versions of checkpoint and devtools
#' install.dates(c('checkpoint', 'devtools'), Sys.Date() - 1)
#'
#' # install yesterday's version of checkpoint and the day before's devtools
#' install.dates(c('checkpoint', 'devtools'), Sys.Date() - 1:2)
#'
#' }
install.dates <- function (pkgs,
                           dates,
                           lib,
                           ...) {

  # number of packages to install
  n_pkgs <- length(pkgs)

  if (!inherits(dates, c('character', 'Date'))) {
    stop ('dates must be a vector of class character or Date')
  }

  if (length(dates) == 1) {
    dates <- rep(dates, n_pkgs)
  }

  if (length(dates) != n_pkgs) {
    stop ('dates must be have either length one, or the same length as pkgs')
  }

  # coerce dates to character
  dates <- as.character(dates)

  # check none of the dates are before the first MRAN date
  if (any(as.Date(dates) < as.Date('2014-09-17'))) {
    stop (sprintf('cannot install packages before 2014-09-17 as this is the
                  earliest date archived on MRAN.
                  Found date: %s',
                  min(as.Date(dates))))
  }

  # loop through packages installing them
  for (i in 1:n_pkgs) {

    # get package and date
    pkg <- pkgs[i]
    date <- dates[i]

    # define repository
    repos <- paste0('https://MRAN.revolutionanalytics.com/snapshot/', date)

    install.packages(pkgs = pkg,
                     lib = lib,
                     repos = repos,
                     ...)

  }

}
