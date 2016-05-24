#' Confirm that all input vectors are the same length.
#'
#' 
#' @param ... Two or more vectors to be compared.
#' @return Returns an error if the vectors are of unequal length, returns a warning if only one vector is supplied, and returns nothing if the vectors are of equal length.
#' @export
#' @examples
#' \dontrun{equal.lengths(rnorm(10), rnorm(10), rnorm(9))}
#' \dontrun{equal.lengths(rnorm(10))}
#' \dontrun{equal.lengths(rnorm(10), rnorm(10), rnorm(10))}
equal.lengths <- function(...){

    vec.list <- list(...)

    # Test that at least one vector was supplied
    if(length(vec.list) < 1){stop("No vectors supplied.")}

    # Warn if only one vector was supplied
    if(length(vec.list) == 1){warning("Only one vector supplied - no basis for comparison of vector lengths.")}

    # Throw error if vectors are of unequal lengths
    if(!all(sapply(vec.list, length) == length(vec.list[[1]]))){
        stop("All input vectors must be the same length.")
    }

}

