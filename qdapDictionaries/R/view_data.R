#' List all data sets available in a \pkg{qdapDictionaries}
#' 
#' Lists and describes all the data sets available in \pkg{qdapDictionaries}.
#' 
#' @param package The name of the package.
#' @seealso \code{\link[utils]{data}}
#' @return Returns the data sets of \pkg{qdapDictionaries} as a dataframe.
#' @export
#' @examples
#' view_data()
view_data <-
function(package = "qdapDictionaries") {
    results <- data(package = package)[["results"]]
    o <- as.data.frame(results[, 3:4], stringsAsFactors = FALSE)
    class(o) <- c("view_data", "data.frame")
    return(o)
}


#' Prints a view_data Object
#' 
#' Prints a view_data object.
#' 
#' @param x The view_data object.
#' @param \ldots ignored
#' @method print view_data
#' @export
print.view_data <-
function(x, ...) {
    width <- options()[["width"]]
	options(width=10000)
    on.exit(options(width=width))
    print(left.just(x))
    return()
}
