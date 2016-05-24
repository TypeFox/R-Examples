## ============================================================================
##
## package - an S3 class that contains attributes about an R package that can
##           be used to load or install that package
##
## ============================================================================

get_repo_info <- function(userRepoSplitSymbol) {
    
    # Returns a function that, given a symbol, breaks a package name and repo
    # combination into its username and repo components if need be
    function(package_name, username_and_repo = FALSE) {
        
        symbol_pos <- get_symbol_pos(package_name, 
                                            userRepoSplitSymbol)
        username <- substr(package_name, 1, (symbol_pos - 1))
        repo <- substr(package_name, (symbol_pos + 1), 
                       nchar(package_name))
        
        list(username = username, repo = repo, 
             userrepo = paste0(username, "/", repo))
    }
}

get_github_info <- get_repo_info("/")
get_bitbucket_info <- get_repo_info("\\$")

## For export -----------------------------------------------------------------

#' package()
#'
#' A package object that contains at minimum the name of a package. If the 
#' package exists on CRAN, the name of the package is used to install that
#' package from CRAN. A forward slash may be used in the package name to 
#' indicate a GitHub username and repo that a package is located in. A dollar
#' sign may be used in the package name to indicate a Bitbucket username and 
#' repo that a package is located in.
#' 
#' If a package is not inteded to be installed from CRAN, Github (public), or 
#' Bitbucket (public) you may optionally supply a function to the 
#' \code{install} argument that installs the package, with additional arguments
#' to the function supplied via the \code{...} argument.
#'
#' @param name the name of a package, or, if the package is on Github, the
#' username and repo of the package, ex. Rdatatable/data.table, where 
#' Rdatatable is the GitHub username, and data.table is the repo name. The 
#' same process works with Bitbucket, except \code{$} is used instead of 
#' \code{/} to separate the username and repo.
#' @param install (optional) a function used to install a package.
#' @param ... (optional) additional arguments to the \code{install} function.
#' Another way of supplying additonal parameters to \code{install} is to use
#' the \code{purrr::partial} function.
#' @param load_type (default is "attach") should the package be loaded, or 
#' attached? See http://r-pkgs.had.co.nz/namespace.html to learn more.
#' @export
#' @examples
#' \dontrun{
#' this_package <- package("dplyr")
#' github_package <- package("Rdatatable/data.table")
#' that_package <- package("jakePackage", devtools::install_bitbucket, 
#' repo = "repo", username = "user", password = "password")
#' }

package <- function(name, install = utils::install.packages, ...,
                    load_type = "attach") {
    
    # Error handling
    arguments <- list(...)
    assert_that(is_function(install))
    if (!load_type %in% c("load", "attach")) 
        stop("load_type may only be 'load' or 'attach'")
    
    # Store the function name (for printing), same for ...
    func_name <- deparse(substitute(install))
    arg_names <- "None"
    if (length(arguments) > 0) 
        arg_names <- "Additional arguments supplied for the install function"
    
    # If :: is appended to the package name, set load_type to "load"
    if (grepl("::", name)) {
        load_type <- "load"
        name <- sub("::", "", name)
    }
    
    # Determine if the user inputted a GitHub package
    if (grepl("/", name)) {
        
        # GitHub
        github_info <- get_github_info(name)
        install <- install_github
        func_name <- "install_github"
        arguments <- list(repo = github_info$userrepo)
        name <- github_info$repo
        
    } else if (grepl("\\$", name)) {
        
        # Bitbucket
        bitbucket_info <- get_bitbucket_info(name)
        install <- install_bitbucket
        func_name <- "install_bitbucket"
        arguments <- list(repo = bitbucket_info$userrepo)
        name <- bitbucket_info$repo
    }
    
    # Create, output the package object
    if (length(arguments) > 0) {
        
        # User supplies install function and separate arguments for that func
        package_obj <- list(name = name, install_func = install,
                            install_args = arguments, load_type = load_type,
                            func_name = func_name, arg_names = arg_names)
    } else if (!all.equal(install, utils::install.packages)) {
        
        # User supplies an install function but no install arguments
        package_obj <- list(name = name, install_func = install,
                            install_args = NULL, load_type = load_type,
                            func_name = func_name, arg_names = arg_names)
    } else {
        
        # Package going to be installed via CRAN
        package_obj <- list(name = name, install_func = "default",
                            install_args = NULL, load_type = load_type,
                            func_name = func_name, arg_names = arg_names)
    }
    
    structure(package_obj, class = "package_obj")
}

#' is.package_obj()
#' Returns TRUE if obj is a package_obj object
#' @param obj an R object
#' @export
is.package_obj <- function(obj) inherits(obj, "package_obj")

#' load_package()
#' Function to load a package given a package object
#' @param pack a \code{package_obj} object 
#' @export
load_package <- function(pack) UseMethod("load_package")

# Loads a package object
#' @export
load_package.package_obj <- function(pack) {
    
    # Either load or attach the package, and return a logical to indicate 
    # whether or not the package loaded/attached correctly
    load_type <- pack$load_type
    if (load_type == "attach") {
        status <- require(pack$name, character.only = TRUE)
        
    } else if (load_type == "load") {
        status <- requireNamespace(pack$name)
        
    } else {
        stop("load_type can only be 'load' or 'attach'")
    }
    
    status
}

#' install_package()
#' Function to install a package given a package object
#' @param pack a \code{package_obj} object 
#' @export
install_package <- function(pack) UseMethod("install_package")

# Installs a package object
#' @export
install_package.package_obj <- function(pack) {
    
    # Construct our install arguments
    if (is.null(pack$install_args)) {
        if (all.equal(pack$install_func, "default")) {
            install_args <- list(pack$name)
            pack$install_func <- utils::install.packages
        } else {
            install_args <- list()
        }
    } else {
        install_args <- list(pack$install_args)
    }
    
    # Try installing the package, but don't throw an error
    x <- try(do.call(pack$install_func, install_args))
             
    # Return TRUE if installation was successful, FALSE otherwise
    if (class(x) == "try-error") {
        return(FALSE)
    } else {
        return(TRUE)
    }
}

# Print method for package_obj
#' @export
print.package_obj <- function(x, ...) {
    
    # Package name
    cat("Package:", x$name, "\n")
    
    # Function name
    if (x$func_name == "utils::install.packages") {
        cat("Install function: (default) utils::install.packages\n")
    } else {
        cat("Install function:", x$func_name, "\n")
    }
    
    # Additional argument names
    if (x$arg_names == "None") {
        cat("Additional install arguments: (default) None\n")
    } else {
        cat("Additional install arguments:", x$arg_names, "\n")
    }
    
    # Load type
    if (x$load_type == "attach") {
        cat("Load/attach: (default) attach\n")
    } else {
        cat("Load/attach:", x$load_type, "\n")
    }
}
