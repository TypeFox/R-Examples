###############################################################################
## one-step estimator
###############################################################################
setMethod("oneStepEstimator", signature(x = "numeric", 
                                        IC = "InfluenceCurve",
                                        start = "numeric"),
    function(x, IC, start){
        nrvalues <- dimension(IC@Curve)
        if(is.list(start)) start <- unlist(start)
        if(nrvalues != length(start))
            stop("length of 'start' != dimension of 'Curve'")
        
        res <- start + rowMeans(evalIC(IC, as.matrix(x)), na.rm = TRUE)

        return(res)        
    })
setMethod("oneStepEstimator", signature(x = "numeric", 
                                        IC = "InfluenceCurve",
                                        start = "list"),
    function(x, IC, start){
        return(oneStepEstimator(x, IC, unlist(start)))
    })
setMethod("oneStepEstimator", signature(x = "matrix", 
                                        IC = "InfluenceCurve",
                                        start = "numeric"),
    function(x, IC, start){
        nrvalues <- dimension(IC@Curve)
        if(is.list(start)) start <- unlist(start)
        if(nrvalues != length(start))
            stop("length of 'start' != dimension of 'Curve'")
        if(ncol(x) != IC@Curve[[1]]@Domain@dimension)
            stop("'x' has wrong dimension")
        
        res <- start + rowMeans(evalIC(IC, x), na.rm = TRUE)

        return(res)        
    })
setMethod("oneStepEstimator", signature(x = "matrix", 
                                        IC = "InfluenceCurve",
                                        start = "list"),
    function(x, IC, start){
        return(oneStepEstimator(x, IC, unlist(start)))
    })
