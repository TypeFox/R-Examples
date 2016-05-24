#' @export
packagesOf <- function(...) UseMethod("packagesOf")

#' Identify the packages of the globals
#'
#' @param globals A Globals object.
#' @param \dots Not used.
#'
#' @return Returns a character vector of package names.
#'
#' @aliases packagesOf
#' @export
packagesOf.Globals <- function(globals, ...) {
  ## Scan 'globals' for which packages needs to be loaded.
  ## This information is in the environment name of the objects.
  pkgs <- sapply(globals, FUN=function(obj) {
    environmentName(environment(obj))
  })

  ## Drop "missing" packages, e.g. globals in globalenv().
  pkgs <- pkgs[nzchar(pkgs)]

  ## Drop global environment
  pkgs <- pkgs[pkgs != "R_GlobalEnv"]

  ## Keep only names matching loaded namespaces
  pkgs <- intersect(pkgs, loadedNamespaces())

  ## Packages to be loaded
  pkgs <- sort(unique(pkgs))

  ## Sanity check
  stopifnot(all(nzchar(pkgs)))

  pkgs
} # packagesOf()

