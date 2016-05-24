#' Creates a package matching the given description and dependencies.
#' @param name the package name
#' @param title the package title
#' @param dependencies a data.frame of package dependencies, including the package names, comparators, and versions
#' @param path the path in which to create the package. Defaults to the current path
#' @return the package, as constructed using the `devtools` `as.package` function
#' @export
#' @examples
#' # Create a simple package with no dependencies:
#' path <- tempdir()
#' name <- 'simplepackage'
#' package <- create_package(name, 'A simple mock package', data.frame(), path)
create_package <- function(name, title, dependencies, path = '.') {

  full_path <- file.path(path, name)

  if(file.exists(full_path)) {
    stop(sprintf("The package directory already exists: [%s]", full_path))
  }

  dir.create(full_path)
  description <- create_package_description(name, title, dependencies)
  write(description, file=file.path(full_path, 'DESCRIPTION'))

  as.package(full_path)

}

