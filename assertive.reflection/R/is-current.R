#' Is this version of R up to date?
#' 
#' Check if this version of R is as new as the current release version of R.
#' @param cran A string giving the URL of the CRAN repository to check.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return An object of class \code{R_system_version} giving the current release
#' version of R.
#' @note Development versions of R can have versions higher than the current
#' release version of R.  For convenience, these will return \code{TRUE}.
#' @examples
#' \donttest{
#' # This example is marked "don't test" since it requires an 
#' # internet connection and is potentially long running
#' is_r_current()
#' }
#' @export
is_r_current <- function(cran = getOption("repos", c(CRAN = "http://cran.r-project.org"))["CRAN"])
{
  this_version <- getRversion()
  current_version <- get_current_r(cran = cran)
  if(this_version < current_version)
  {
    return(
      false(
        gettext("This version of R is %s but the current version is %s."), 
        this_version,
        current_version
      )
    )
  }
  TRUE
}

is_current_r <- function(cran = getOption("repos", c(CRAN = "http://cran.r-project.org"))["CRAN"])
{
  .Deprecated("is_r_current")
  is_r_current(cran = cran)
}

#' Is the installed version of a package current?
#' 
#' Checks to see if the installed version of a package is current.
#' @param x A string giving a package name.
#' @param lib.loc A character vector of paths to local package libraries.
#' @param repos A character vector of URLs to repositories to check for new
#' package versions.
#' @param type Check the repository for source or binary packages?
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @param .xname Not intended to be used directly.
#' @return The \code{is_*} functions return \code{TRUE} or \code{FALSE}.
#' The \code{assert_*} functions throw an error in the event of failure.
#' @seealso \code{\link[utils]{old.packages}}, on which this is based, which
#' has advanced usage features.
#' @examples 
#' \donttest{
#' # This test is marked "dont-test" since it involves a connection to 
#' # repositories which is potentially long running.
#' is_package_current("assertive.reflection")
#' }
#' @importFrom utils installed.packages
#' @importFrom utils old.packages
#' @importFrom assertive.base coerce_to
#' @importFrom assertive.base use_first
#' @export
is_package_current <- function(x, lib.loc = .libPaths(), 
  repos = getOption("repos"), type = getOption("pkgType"), 
  .xname = get_name_in_parent(x))
{
  # TODO: what is the behaviour when the package is installed with 
  # different versions in multiple local libraries?
  
  x <- coerce_to(use_first(x), "character", .xname)
  ip <- installed.packages(lib.loc = lib.loc)
  requestedPkgIsInstalled <- x %in% rownames(ip)
  if(!all(requestedPkgIsInstalled))
  {
    stop(
      "The following packages are not installed and cannot be checked: ", 
      toString(x[!requestedPkgIsInstalled])
    )
  }
  # Prevent "trying to use CRAN without setting a mirror" in contrib.url
  if ("@CRAN@" %in% repos)
  {
    # It should only be the CRAN repo assigned the dummy value "@CRAN@"
    # If users have assigned it to other repos, then this next line is silly,
    # but that should be so rare as to not need worrying about.
    repos[repos == "@CRAN@"] <- "https://cloud.r-project.org"
  }
  op <- old.packages(
    lib.loc = lib.loc,
    repos = repos,
    type = type,
    instPkgs = ip[x, , drop = FALSE])
  if(!is.null(op))
  {
    return(
      false(
        gettext("%s is out of date.  The installed version is %s but the latest %s version is %s."),
        x,
        as.character(op[, "Installed"]),
        type,
        as.character(op[, "ReposVer"])
      )
    )
  }
  TRUE
}


# The smart implementation of this function uses rvest, but we don't want 
# the dependency, so use readLines + regex matching instead.
# get_current_r <- function(cran = getOption("repos", c(CRAN = "http://cran.r-project.org"))["CRAN"])
# {
#   doc <- rvest::html(paste(cran, "sources.html", sep = "/"))
#   # Version should be contained in the first link on this page
#   `%>%` <- rvest::`%>%`
#   version_string <- doc %>% 
#     rvest::html_node("li > a") %>%
#     rvest::html_text()
#   R_system_version(substring(version_string, 3, nchar(version_string) - 7))
# }

get_current_r <- function(cran = getOption("repos", c(CRAN = "http://cran.r-project.org"))["CRAN"])
{
  if(cran == "@CRAN@")
  {
    cran <- "http://cran.r-project.org"
  }
  lines <- readLines(paste(cran, "sources.html", sep = "/"))
  rx <- "R-(\\d\\.\\d\\.\\d)"
  version_string <- regmatches(lines, regexpr(rx, lines))
  R_system_version(substring(version_string, 3))
}

