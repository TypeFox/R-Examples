#' packages()
#'
#' Loads one or more packages, and installs them after a user prompt if they
#' are not already installed. You may supply either package names or 
#' \code{package_obj} objects to this function. You can install Github packages
#' by supplying \code{username/repo} to this function, or \code{username$repo}
#' for Bitbucket packages.
#' 
#' \code{package_obj} allows you to supply it an install function if
#' the package isn't on CRAN or in a public GitHub or Bitbucket repo.
#' 
#' You may also add \code{::} to the end of a package name to load
#' the package instead of attaching the package (this means that the package
#' will not be added to the search list, but will be available to access via
#' the \code{::} operator).
#' 
#' @param ... one or more package names or \code{package_obj} objects
#' @param prompt (default is TRUE) prompt the user before installing packages?
#' For interactive use keeping the default makes the most sense. For 
#' interactive scripts, or scripts that you are sharing with those you trust,
#' it may make more sense to switch this to FALSE.
#' @export
#' @examples
#' \dontrun{
#' packages("dplyr", "ggplot2", "rvest", "magrittr")
#' packages("dplyr::", "Rdatatable/data.table")
#' }

packages <- function(..., prompt = TRUE) {
    
    # Gather ..., turn any non-package_obj objects into package_obj objects
    arguments <- list(...)
    packages <- create_package_objs(arguments)
    
    # Determine which packages need to be installed
    install_these <- list()
    for (i in seq_along(packages))
        if (!is_package_installed(packages[[i]]$name))
            install_these[[packages[[i]]$name]] <- packages[[i]]
    
    # Prompt the user to install uninstalled packages
    if ((length(install_these) > 0) & prompt) {
        
        if (length(install_these) == 1) {
            
            # One package needs to be installed
            message("The following package is not installed: ", 
                    paste(names(install_these), collapse = ", "))
            message("Would you like to install this package now? [y/n]")
            command <- scan(what = character(), n = 1, quiet = TRUE)
            if (command != "y")
                stop("You choose not to install this package. Please choose y", 
                     " or install this package yourself to proceed.")
        } else {
            
            # Two or more packages need to be installed
            message("The following packages are not installed: ", 
                    paste(names(install_these), collapse = ", "))
            message("Would you like to install these packages now? [y/n]")
            command <- scan(what = character(), n = 1, quiet = TRUE)
            if (command != "y")
                stop("You choose not to install these packages. Please choose",
                     " y or install these packages yourself to proceed.")
        }
        
        # Install the packages that are not installed
        do.call(install_packages, install_these)
        
    } else if (length(install_these) > 0) {
    
        # Install the packages that are not installed
        do.call(install_packages, install_these)
    }
    
    # Load all of the packages
    do.call(libraries, packages)
}
