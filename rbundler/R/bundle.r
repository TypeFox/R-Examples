#' Bundles a package and it's dependencies into a library.
#'
#' Dependencies are installed into the package's bundle library. The
#' library is also added to this session's .libPaths.
#'
#' Note that repository and pkgType options are temporarily overridden,
#' according to the user's options, and set back to their previous values after
#' bundle completes.
#' @param pkg package description, can be path or package name.
#' @param bundle_path path to the bundle. Defaults to '.Rbundle' under the package directory
#' @param overwrite whether to delete the existing bundle library and re-install all packages. This can be necessary when
#'        upgrading or downgrading package dependencies. Defaults to FALSE
#' @param dependencies which package dependencies to install. Defaults to c("Depends", "Imports", "LinkingTo", "Suggests")
#' @importFrom devtools install
#' @importFrom devtools as.package
#' @export
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
bundle <- function(pkg = '.', bundle_path = file.path(pkg, '.Rbundle'), overwrite = FALSE,
  dependencies = c("Depends", "Imports", "LinkingTo", "Suggests")) {

  package <- as.package(pkg)

  if(overwrite && file.exists(bundle_path)) {
    message(sprintf("Overwriting existing bundle at [%s]", bundle_path))
    unlink(bundle_path, recursive=TRUE)
  }

  dir.create(bundle_path, recursive = TRUE, showWarnings = FALSE)

  r_libs_user <- construct_r_libs_user(bundle_path)

  update_renviron_file(path = package$path, r_libs_user = r_libs_user)
  update_current_environment(lib = bundle_path, r_libs_user = r_libs_user)

  message("Bundling package ", package$path, " dependencies into library ", bundle_path)

  if (!is.null(package$depends)) {

    depends <- devtools::parse_deps(package$depends)

    if(nrow(depends) > 0) {
      apply(
        depends,
        1,
        FUN = function(d) {
          install_version(
            d['name'], d['version'], d['compare'],
            dependencies = dependencies
          )
        }
      )
    }

  }

  install(pkg)

  invisible()

}

#' Constructs a new R_LIBS_USER setting using the current libraries and the new bundle library.
#' @param bundle_path the new bundle path
#' @return r_libs_user colon-separated libraries
construct_r_libs_user <- function(bundle_path) {
  current_libs <- unlist(strsplit(Sys.getenv('R_LIBS_USER'), split=':'))
  new_libs <- unique(c(current_libs, bundle_path))
  paste(new_libs, collapse=':')
}

#' Updates the current environment.
#' @export
#' @param lib the R library to add.
#' @param r_libs_user the new value of R_LIBS_USER
update_current_environment <- function(lib, r_libs_user) {

  Sys.setenv(R_LIBS_USER=r_libs_user)
  .libPaths(c(lib, .libPaths()))

  message("R_LIBS_USER=", r_libs_user)
  message(".libPaths() =", .libPaths())

  invisible()

}

#' Updates a .Renviron file in the given path.
#' @export
#' @param path to the .Renviron file
#' @param r_libs_user the new value of R_LIBS_USER
update_renviron_file <- function(path, r_libs_user) {

  renviron <- file(file.path(path, ".Renviron"))
  writeLines(sprintf("R_LIBS_USER='%s'", r_libs_user), renviron)
  close(renviron)

  invisible()

}

