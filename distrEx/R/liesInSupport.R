setMethod("liesInSupport", signature(object = "DiscreteMVDistribution",
                                     x = "numeric"),
    function(object, x){
        k <- dimension(img(object))
        if(length(x) != k)
            stop("'x' has wrong dimension") 
        supp <- support(object)
        ind <- colSums(apply(supp, 1, "==", x)) == k

        return(any(ind))
    })
setMethod("liesInSupport", signature(object = "DiscreteMVDistribution",
                                     x = "matrix"),
    function(object, x){ 
        if(ncol(x) != dimension(img(object)))
            stop("'x' has wrong dimension") 

        res <- logical(nrow(x))
        for(i in 1:nrow(x)) res[i] <- liesInSupport(object, x[i,])

        return(res)
    })
