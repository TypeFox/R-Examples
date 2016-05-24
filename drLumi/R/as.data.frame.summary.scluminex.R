#' @export
as.data.frame.summary.scluminex <- function(x,row.names = NULL, 
                                            optional = FALSE, ...) {
    if (inherits(x,"summary.scluminex")){
        return(data.frame(unclass(x),...))
    } else {
        stop("'x' is not an object of class summary.scluminex")
    }
}
