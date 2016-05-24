#' libraries()
#'
#' A vectorized version of the \code{library} function that accepts one or
#' more package names to call \code{library} on. Unlike \code{library},
#' the \code{libraries} function does not use non-standard evaluation (NSE)
#' on its inputs, meaning that you must supply your package names as strings,
#' or as variables.
#' 
#' You may supply either package names and \code{package_obj} objects to this
#' function. You may also add \code{::} to the end of a package name to load
#' the package instead of attaching the package (this means that the package
#' will not be added to the search list, but will be available to access via
#' the \code{::} operator).
#' 
#' @param ... one or more package names or \code{package_obj} objects
#' @export
#' @examples
#' \dontrun{
#' libraries("dplyr", "ggplot2", "rvest", "magrittr")
#' libraries("dplyr::", "Rdatatable/data.table")
#' }

libraries <- function(...) {

    # Gather ..., turn any non-package_obj objects into package_obj objects
    arguments <- list(...)
    packages <- create_package_objs(arguments)
    
    # Load each package
    packages_loaded <- vapply(packages, load_package, logical(1))
    
    # Return whether or not the packages all loaded successfully
    if (all(packages_loaded)) {
        message("All packages loaded successfully")
    } else {
        message("\n--------------------------------------------------")
        message("The following packages did not load successfully: ")
        message(names(packages_loaded[!packages_loaded]), "\n")
    }
}
