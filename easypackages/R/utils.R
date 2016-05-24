## ============================================================================
##
## Utilities for the packages package, including assertations
##
## ============================================================================

## Utility functions ----------------------------------------------------------

get_symbol_pos <- function(string, symbol) {
    
    # Returns the position of a symbol in a string
    gregexpr(symbol, string)[[1]][1]
}

se_lapply <- function(X, FUN, ...) {
    
    # "Side effect" only version of lapply - nothing is returned, the function
    # is simply called, and therefore produces a side effect
    for (i in seq_along(X)) {
        FUN(X[[i]], ...)
    }
}

# Names of all installed packages
installed_packages <- function() unname(utils::installed.packages()[, 1])

# Returns TRUE if a package is installed, FALSE otherwise
is_package_installed <- function(package_name)
    package_name %in% installed_packages()

create_package_obj <- function(argument) {
    
    # Create a single package object if it isn't already a package object
    if (is.package_obj(argument)) {
        mypackage <- argument
    } else {
        mypackage <- package(argument)
    }
    
    mypackage
}

lowest_level <- function(x) {
    
    # Is an individual element at the lowest level? Lowest level means that
    # we can call create_all_package_objs() on an object
    if (is.package_obj(x) | (length(x) == 1 & is.character(x))) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

create_package_objs <- function(mylist) {
    
    # Go down the list recursively, turn non-package_objs into package_objs
    packagez <- list()
    for (i in seq_along(mylist)) {
        if (lowest_level(mylist[[i]])) {
            mypackage <- create_package_obj(mylist[[i]])
            packagez <- c(packagez, list(mypackage))
        } else {
            packagez <- c(packagez, create_package_objs (mylist[[i]]))
        }
    }
    
    packagez
}

## Import functions -----------------------------------------------------------

#' @importFrom assertthat assert_that
#' @importFrom assertthat is.flag
#' @importFrom devtools install_github
#' @importFrom devtools install_bitbucket

## Assertations ---------------------------------------------------------------

package_installed <- function(package) is_package_installed(package)
assertthat::on_failure(package_installed) <- function(call, env) {
    paste0("Package '", deparse(call$package), "' is not installed.")
}

is_function <- function(x) is.function(x)
assertthat::on_failure(is_function) <- function(call, env) {
    paste0("Argument '", deparse(call$x), "' must be a function, but isn't.")
}

is_vector <- function(x) 
    mode(x) %in% c("logical", "numeric", "complex", "character")

assertthat::on_failure(is_vector) <- function(call, env) {
    paste0("Argument '", deparse(call$x), "' is not a vector.")
}

no_duplicates <- function(x) sum(duplicated(x)) == 0
assertthat::on_failure(no_duplicates) <- function(call, env)
    paste0("Duplicates exist in ", deparse(call$x), " when they should be none.")
