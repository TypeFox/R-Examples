#' listWithDefaults
#' 
#' Takes arguments as \code{base::list} to create a list.
#' If any arguments are present in \code{defaultList}, but absent in ..., then the values from \code{defaultList} are used.
#' 
#' @note Argument order is not controlled.  Non-default arguments come first in the order specified followed by all default arguments.
#' @export
#' @importFrom assertthat assert_that
#' @param ... objects, \emph{must} be named
#' @param defaultList a \emph{named} list containing the default values
#' @examples 
#' listWithDefaults(defaultList=list(a=2,b=2))
#' listWithDefaults(a=1,defaultList=list(a=2,b=2))
#' listWithDefaults(b=1,defaultList=list(a=2,b=2))
#' listWithDefaults(a=1,b=1,defaultList=list(a=2,b=2))

listWithDefaults <-
function(..., defaultList=NULL) {
    result <- list(...)
    assert_that(length(unique(names(result)))==length(names(result)))
    assert_that(length(unique(names(defaultList)))==length(names(defaultList)))
    assert_that(!is.null(defaultList))
    assert_that(all(names(result) %in% names(defaultList)))
    for (thisName in names(defaultList)[!names(defaultList) %in% names(result)]) {
        result[thisName] <- defaultList[thisName] 
    }
    return(result)
}

