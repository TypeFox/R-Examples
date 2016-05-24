setMethod("getDistance", "CovNA", function(obj){
    if(!is(obj@mah,"NULL"))
        return(obj@mah)

    if(is(getData(obj), "NULL"))
        stop("Cannot compute distances: no data provided")

    ## dd <- mahalanobis(obj@X, obj@center, obj@cov)
    dd <- .mah.na(obj@X, getCenter(obj), getCov(obj))

    ## check if the method is called on an object - i.e. cc=CovMcd(xxx); getDistance(cc)
    ##  or on the constructor call - i.e. getDistance(CovMcd(xxx))
    ## Do nothing in the second case
    ## Our expression is 'obj@mah <- dd' and the parse tree is
    ##  (as a string array) will be c("<-", "obj@mah", "dd")
    ## We are interested if there are () in the second element
    ##
    expr <- substitute(obj@mah <- dd)
    ss <- as.character(expr)
    if(length(grep(pattern="(", x=ss[2], fixed=TRUE)) == 0)
        eval.parent(expr)

    return(dd)
})

setMethod("getFlag", "CovNA", function(obj, prob=0.975){
    if(!is(obj@flag,"NULL") && missing(prob))
        return(obj@flag)

    p <- ncol(getCov(obj))

    if(!is(obj@mah,"NULL"))
        dd <- obj@mah
    else if(is(getData(obj), "NULL"))
        stop("Cannot compute distances: no data provided")
    else
        dd <- .mah.na(obj@X, getCenter(obj), getCov(obj))

    fl <- .dflag(dd$d, prob, dd$pp)

    ## check if the method is called on an object - i.e. cc=CovMcd(xxx); getFlag(cc)
    ##  or on the constructor call - i.e. getFlag(CovMcd(xxx))
    ## Do nothing in the second case
    ## Our expression is 'obj@flag <- fl' and the parse tree is
    ##  (as a string array) will be c("<-", "obj@flag", "fl")
    ## We are interested if there are () in the second element
    ##
    expr <- substitute(obj@flag <- fl)
    ss <- as.character(expr)
    if(length(grep(pattern="(", x=ss[2], fixed=TRUE)) == 0)
        eval.parent(expr)
    return(fl)
})



setMethod("summary", "CovNA", function(object, ...){

    new("SummaryCovNA", covobj=object, evals=eigen(object@cov)$values)

})

setMethod("show", "SummaryCovNA", function(object){

    cat("\nCall:\n")
    print(object@covobj@call)

    digits = max(3, getOption("digits") - 3)
    cat("\nEstimate of Location: \n")
    print.default(format(getCenter(object), digits = digits), print.gap = 2, quote = FALSE)
    cat("\nEstimate of Covariance: \n")
    print.default(format(getCov(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nEigenvalues of covariance matrix: \n")
    print.default(format(getEvals(object), digits = digits), print.gap = 2, quote = FALSE)

    cat("\nMahalanobis Distances: \n")
    dd <- getDistance(object)
    print.default(format(as.vector(dd$d), digits = digits), print.gap = 2, quote = FALSE)
})

##setMethod("impute", "CovNA", function(obj, ...){
##    impWins(obj@X, obj@center, obj@cov, !getFlag(obj))
##})
##setMethod("impute", "CovNAMcd", function(obj, ...){
##    impWins(obj@X, obj@center, obj@cov, !getFlag(obj))
##})
