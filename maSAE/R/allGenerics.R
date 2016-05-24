#' @export
#' @docType methods
#' @rdname predict-methods
#' @param object a model object for which prediction is desired. 
#' @param ... Arguments to be passed to methods.
setGeneric(name = "predict"
            , def = function(object, ...){standardGeneric("predict")}
            )

