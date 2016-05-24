
# dataSet####
setAs("dataSet", "matrix", function(from) as(from@data, "matrix"))

setMethod("dim", signature(x = "dataSet"), function(x) dim(x@data))

setMethod("ncol", signature(x = "dataSet"), function(x) ncol(x@data))

setMethod("nrow", signature(x = "dataSet"), function(x) nrow(x@data))

setMethod("colRatings", signature(x = "dataSet"), function(x) {
    s <- apply(x@data, 2, function(i) sum(i != 0))
    names(s) <- colnames(x@data)
    s
})


setMethod("rowRatings", signature(x = "dataSet"), function(x) {
    s <- apply(x@data, 1, function(i) sum(i != 0))
    names(s) <- rownames(x@data)
    s
})

setMethod("numRatings", signature(x = "dataSet"), function(x) {
    sum(x@data != 0)
})


setMethod("[", signature(x = "dataSet", i = "ANY", j = "ANY", drop = "missing"), function(x, i, j, ..., drop) {
    
    if (missing(i)) 
        i <- 1:nrow(x)
    if (missing(j)) 
        j <- 1:ncol(x)
    
    x@data <- x@data[i, j, ..., drop = FALSE]
    x
})
setMethod("sparsity", signature(x = "dataSet"), function(x) {
    
    spars <- 1 - numRatings(x)/(nrow(x@data) * ncol(x@data))
    spars
})

setMethod("show", signature(object = "evalModel"), function(object) {
    if (object@data@binary) {
        cat(object@folds, "- fold cross validation model on the binary dataset with", nrow(object@data), "users and", ncol(object@data), "items.")
    } else {
        cat(object@folds, "- fold cross validation model on the dataset with", nrow(object@data), "users and", ncol(object@data), "items.")
        
    }
})


# recResultsClass####
setMethod("show", signature(object = "recResultsClass"), function(object) {
    print(object@recommended)
})

setMethod("[", signature(x = "recResultsClass", i = "ANY", j = "missing", drop = "missing"), function(x, i, j, ..., drop) {
    
    if (missing(i)) 
        i <- 1:nrow(x)
    
    
    x@recommended[i]
})




# show methods for the trained model####
setMethod("show", signature(object = "SVDclass"), function(object) {
    cat("The model was trained on the dataset using ", object@alg, "algorithm.\nThe algorithm was configured with the following parameters:\n")
    print(as.data.frame(object@parameters))
})

setMethod("show", signature(object = "IBclass"), function(object) {
    cat("The model was trained on the dataset using ", object@alg, "algorithm. \nThe algorithm was configured with the following neighborhood width:", object@neigh)
})

setMethod("show", signature(object = "wALSclass"), function(object) {
    cat("The model was trained on the dataset using ", object@alg, "algorithm. \nThe algorithm was configured with the following parameters:\n")
    print(as.data.frame(object@parameters))
})

setMethod("show", signature(object = "BPRclass"), function(object) {
    cat("The model was trained on the dataset using ", object@alg, "algorithm. \nThe algorithm was configured with the following parameters:\n")
    print(as.data.frame(object@parameters))
})

setMethod("show", signature(object = "PPLclass"), function(object) {
    cat("The model was trained on the dataset using ", object@alg, "algorithm.")
})


setMethod("show", signature(object = "algAverageClass"), function(object) {
    cat("The model was trained on the dataset using ", object@alg, "algorithm.")
}) 
