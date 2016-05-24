#' Rbundler is an R package dependency management utility.
#'
#' @name rbundler
#' @title A package dependency management utility.
#' @docType package
#' @author Yoni Ben-Meshulam \email{yoni.bmesh@@gmail.com}
#' @examples
#'\dontrun{
#' # Run bundle in the current path:
#' bundle()
#' # Check for the new `.Rbundle` entry in `.libPaths`:
#' .libPaths()
#'
#' lib <- file.path(tempdir(), 'my_bundle_lib')
#' # Run bundle in the current path, overriding the target library:
#' bundle('.', lib)
#' # Check for the new entry in `.libPaths`:
#' .libPaths()
#'}
NULL
