setMethod("rrecsys", signature = c("dataSet"), function(data, alg, ...) {
    recom <- rrecsysRegistry$get_entry(alg = alg)
    if (is.null(recom)) 
        stop(paste("Wrong method!!!"))
    
    recom$fun(data, ...)
}) 
