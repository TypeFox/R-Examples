setMethod("[", signature(x="APResult", i="index", j="missing", drop="missing"),
    function(x, i, j, drop=FALSE)
    {
        x@clusters[i]
    })

setMethod("[[", signature(x="APResult", i="index", j="missing"),
    function(x, i, j)
    {
        x@clusters[[i]]
    })

setMethod("[", signature(x="ExClust", i="index", j="missing", drop="missing"),
    function(x, i, j, drop=FALSE)
    {
        x@clusters[i]
    })

setMethod("[[", signature(x="ExClust", i="index", j="missing"),
    function(x, i, j)
    {
        x@clusters[[i]]
    })

setMethod("[", signature(x="AggExResult", i="index", j="missing",
                         drop="missing"),
    function(x, i, j, drop=FALSE)
    {
        lapply(i, function(index) cutree(x, k=index))
    })

setMethod("[[", signature(x="AggExResult", i="index", j="missing"),
    function(x, i, j)
    {
        cutree(x, k=i)
    })

setMethod("similarity", signature(x="APResult"), function(x) x@sim)

setMethod("similarity", signature(x="AggExResult"), function(x) x@sim)

setMethod("similarity", signature(x="ExClust"), function(x) x@sim)
