#' Check, Install and Load CRAN & USGS GRAN Packages
#'
#' Checks the local R library(ies) to see if the required package(s)
#'	is/are installed or not. If the package(s) is/are not installed,
#'	then the package(s) will be installed along with the required
#'	dependency(ies). This function pulls source or binary packages
#'	from the Revolution Analytics CRAN mirror and/or the USGS GRAN
#'	Repository. Lastly, the chosen package(s) is/are loaded.
#'
#' @param package1 R package.
#' @param ... Any additional R packages.
#'
#' @return Install and load \code{package1} and other packages listed
#'	in the arguments \code{...} and all of the dependencies.
#'
#'
#' @author maloneypatr, Irucka Embry
#'
#' @source Check if R package is installed then load library answered by
#'	maloneypatr is the source for the original function. See
#'	\url{http://stackoverflow.com/questions/15155814/check-if-r-package-is-installed-then-load-library}.
#'
#'
#'
#' @seealso \code{\link[needs]{needs}}, \code{\link[easypackages]{install_packages}},
#'   \code{\link[utils]{installed.packages}}, \code{\link[utils]{install.packages}}
#'
#'
#'
#' @examples
#' \dontrun{
#' library(install.load)
#' # This will install, if not already installed, and load the package(s)
#' install_load("chron")
#' install_load("chron", "data.table", "dataRetrieval", "EGRET", "dplyr")
#' }
#'
#'
#' @import utils
#'
#' @export
# Source begins
install_load <- function (package1, ...) {

  # convert arguments to vector
  packages <- c(package1, ...)

  # start loop to determine if each package is installed
  for (package in packages) {

    # if package is installed locally, load
  if (package %in% rownames(installed.packages()))
      do.call("library", list(package))

    # if package is not installed locally, download and then load
    else {

      install.packages(package, repos =
        c("https://cran.revolutionanalytics.com/", "http://owi.usgs.gov/R/"),
        dependencies = NA, type = getOption("pkgType"))
      do.call("library", list(package))
    }
  }
}
# Source ends
