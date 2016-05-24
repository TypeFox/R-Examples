#' install.versions
#'
#' @description Download and install named versions of packages hosted on
#'  CRAN from the MRAN server.
#'
#' @param pkgs character vector of the names of packages that should be
#'  downloaded and installed
#' @param versions character vector of the versions of packages to be
#'  downloaded and installed. If this has the same length as \code{pkgs}
#'  versions will correspond to those packages. If this has length one
#'  the same version will be used for all packages. If it has any other
#'  length an error will be thrown.
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
#' @name install.versions
#'
#' @examples
#'
#'\dontrun{
#'
#' # install an earlier version of checkpoint
#' install.versions('checkpoint', '0.3.3')
#'
#' # install earlier versions of checkpoint and devtools
#' install.versions(c('checkpoint', 'devtools'), c('0.3.3', '1.6.1'))
#'
#'}
install.versions <- function (pkgs,
                              versions,
                              lib,
                              ...) {

  # number of packages to install
  n_pkgs <- length(pkgs)

  if (!inherits(versions, 'character')) {
    stop ('versions must be a character vector')
  }

  if (length(versions) == 1) {
    versions <- rep(versions, n_pkgs)
  }

  if (length(versions) != n_pkgs) {
    stop ('versions must be have either length one, or the same length as pkgs')
  }

  # get date corresponding to version for this package
  date <- version2date(pkgs = pkgs,
                       versions = versions)

  # install by date
  install.dates(pkgs = pkgs,
                dates = date,
                lib = lib,
                ...)

}
