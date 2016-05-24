#' A wrapper for \code{library}.
#'
#' Tries to load packages. If the packages are not found, they will be installed from
#' the default repository. This function is intended for use in interactive sessions
#' and should not be used by other packages.
#'
#' @param ... [any]\cr
#'   Package names.
#' @return [\code{logical}]: Named logical vector determining the success
#'   of package load.
#' @export
#' @examples \dontrun{
#' lib("BBmisc", "MASS", "rpart")
#' }
lib = function(...) {
  getLib = function(pkg) {
    ok = suppressWarnings(require(pkg, character.only = TRUE))
    if (!ok && !is.error(try(install.packages(pkg)))) {
      ok = require(pkg, character.only = TRUE)
    }
    ok
  }

  pkgs = unique(c(...))
  assertCharacter(pkgs, any.missing = FALSE)
  vlapply(pkgs, getLib)
}
