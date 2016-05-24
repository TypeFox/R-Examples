#' @title versions: Query and Install Specific Versions of Packages on CRAN
#'
#' @name versions-package
#' @description Installs specified versions of R packages
#' hosted on CRAN and provides functions to list available versions and the
#' versions of currently installed packages. These tools can be used to help
#' make R projects and packages more reproducible.
#' \code{versions} fits in the narrow gap between the devtools
#' \code{install_version} function and the \code{checkpoint} package.
#'
#' \code{devtools::install_version} installs a stated package version from
#' source files stored on the CRAN archives. However CRAN does not store
#' binary versions of packages so Windows users need to have RTools installed
#' and Windows and OSX users get longer installation times.
#'
#' \code{checkpoint} uses the Revolution Analytics MRAN server to
#' install packages (from source or binary) as they were available on
#' a given date. It also provides a helpful interface to detect the packages
#' in use in a directory and install all of those packages for a given date.
#' \code{checkpoint} doesn't provide \code{install.packages}-like functionality
#' however, and that's what \code{versions} aims to do, by querying MRAN.
#'
#' As MRAN only goes back to 2014-09-17, \code{versions} can't install packages
#' from before this date.
#'
#' The available functions are:
#' \itemize{
#'  \item \code{\link{available.versions}}
#'  \item \code{\link{install.versions}}
#'  \item \code{\link{install.dates}}
#'  \item \code{\link{installed.versions}}
#' }
#'
#' @docType package
#'
#' @importFrom utils install.packages download.file
#'
#' @examples
#'
#' \dontrun{
#'
#' # list the available versions of checkpoint
#' available.versions('checkpoint')
#'
#' # install a specific version
#' install.versions('checkpoint', '0.3.9')
#'
#' # check the installed version
#' installed.versions('versions')
#'
#' # install checkpoint as of a specific date
#' install.dates('checkpoint', '2014-12-25')
#'
#' }
#'
#'

NULL
