"fitted.grouped" <-
function(object, ...){
    if(!inherits(object, "grouped"))
        stop("Use only with 'grouped' objects.\n")
    object$fitted
}

