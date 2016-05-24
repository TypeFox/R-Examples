#' from_import()
#'
#' Imports one or more functions from a package into the parent environment,
#' overwriting any existing objects in the parent environment sharing the name
#' of the function(s). The package is loaded but not attached to the search
#' list. Inspired by Python.
#' 
#' @param package the package name (length one character vector)
#' @param ... one or more function names (as-is, or as a list)
#' @param verbose give a warning when overwriting existing objects? TRUE by
#' default. 
#' @export
#' @examples
#' \dontrun{
#' from_import("dplyr", "select", "filter")
#' from_import("dplyr", list("select", "filter"))
#' }

from_import <- function(package, ..., verbose = TRUE) {
    
    # Get a list of functions to import
    functions <- list(...)
    if (length(functions) == 1 && is.list(functions[[1]])) {
        functions <- functions[[1]]
    }
    
    # Error handling
    assert_that(package_installed(package))
    assert_that(length(functions) > 0)
    
    # Loop over the functions, assign the function to the global environment
    for (func_name in functions) {
        
        # If the function name exists in the parent environment
        if (verbose) {
            if (exists(func_name, envir = sys.frame(-1))) {
                warning("Object: ", func_name, " already exists in the ",
                        "parent environment. It was just overwritten.")
            }
        }
        
        # Assign the function
        func <- NULL
        eval(parse(text = paste0("func <- ", package, "::", func_name)))
        assign(func_name, func, envir = sys.frame(-1))
    }
}
