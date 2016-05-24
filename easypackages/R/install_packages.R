#' install_packages()
#'
#' Installs one or more packages. Similar to \code{utils::install.packages}, but
#' you may supply either package names and \code{package_obj} objects to this
#' function. You can install Github packages by supplying \code{username/repo}
#' to this function, or \code{username$repo} for Bitbucket packages.
#' 
#' \code{package_obj} allows you to supply it an install function if
#' the package isn't on CRAN or in a public GitHub or Bitbucket repo.
#' 
#' @param ... one or more package names or \code{package_obj} objects
#' @export
#' @examples
#' \dontrun{
#' install_packages("dplyr", "ggplot2", "rvest", "magrittr")
#' }

install_packages <- function(...) {
    
    # Gather ..., turn any non-package_obj objects into package_obj objects
    arguments <- list(...)
    packages <- create_package_objs(arguments)
    
    # Install each package
    packages_installed <- vapply(packages, install_package, logical(1))
    
    # Return whether or not the packages all loaded successfully
    if (all(packages_installed)) {
        message("All packages installed successfully")
    } else {
        message("\n--------------------------------------------------")
        message("The following packages did not install successfully: ")
        message(names(packages_installed[!packages_installed]), "\n")
    }
}
