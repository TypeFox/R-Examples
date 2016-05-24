overlay <- function(modelList = NULL, showLegend = TRUE, Cols = NULL, ...){
    
    if(is.null(modelList))
        stop("You just provided an empty list.")
    isValid <- sapply(modelList, function(tmp) inherits(tmp, "nplr"))

    if(!all(isValid))
        stop(sprintf("model #%s is not an instance of class 'nplr'", which(!isValid)))
    
    .multiCurve(modelList, showLegend, Cols, ...)
}
