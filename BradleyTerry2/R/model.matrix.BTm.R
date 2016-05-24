model.matrix.BTm <- function(object, ...){
    ## set contrasts to NULL as apply to player formula not dummy formula
    object$contrasts <- NULL
    NextMethod("model.matrix")
}
