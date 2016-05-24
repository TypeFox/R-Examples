#' @export
as.data.frame.scluminex <- function(x, row.names = NULL, 
                                    optional = FALSE, ...){
    data <- ldply(lapply(x,function(y){y$data}),rbind)
    data$.id <- NULL
    return(data)
}
