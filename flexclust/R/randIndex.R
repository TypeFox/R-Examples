setGeneric("randIndex", function(x, y, correct=TRUE, original=!correct)
           standardGeneric("randIndex"))

setMethod("randIndex", signature(x="ANY", y="ANY"),
function(x, y, correct=TRUE, original=!correct){
    if(correct)
        comPart(x, y, type="ARI")
    else
        comPart(x, y, type="RI")
})

setMethod("randIndex", signature(x="table", y="missing"),
doRandIndex <- function(x, y, correct=TRUE, original=!correct)
{
    if(length(dim(x))!=2)
        stop("Argument x needs to be a 2-dimensional table.")
    
    n <- sum(x)
    ni <- apply(x, 1, sum)
    nj <- apply(x, 2, sum)
    n2 <- choose(n, 2)

    rand <- NULL
    if(correct){
        nis2 <- sum(choose(ni[ni > 1], 2))
        njs2 <- sum(choose(nj[nj > 1], 2))
        rand <- c(ARI=c(sum(choose(x[x > 1], 2)) -
                  (nis2 * njs2)/n2)/((nis2 + njs2)/2 - (nis2 * njs2)/n2))
    }

    if(original){
        rand <- c(rand, RI = 1 + (sum(x^2) - (sum(ni^2) + sum(nj^2))/2)/n2)
    }

    return(rand)
})

###**********************************************************

setGeneric("comPart", function(x, y, type=c("ARI","RI","J","FM"))
    standardGeneric("comPart"))

setMethod("comPart", signature(x="flexclust", y="flexclust"),
function(x, y, type){
    doComPart(clusters(x), clusters(y), type)
})

setMethod("comPart", signature(x="flexclust", y="numeric"),
function(x, y, type){
    doComPart(clusters(x), y, type)
})

setMethod("comPart", signature(x="numeric", y="flexclust"),
function(x, y, type){
    doComPart(x, clusters(y), type)
})

setMethod("comPart", signature(x="numeric", y="numeric"),
doComPart <- function(x, y, type=c("ARI","RI","J","FM"))
{
    type <- toupper(type)
    if(length(x)!=length(y))
        stop("x an y must have the same length")

    nxx <- countPairs(x, y)

    res <- NULL
    if("ARI" %in% type)
        res <- c(doRandIndex(table(x,y), correct=TRUE))
    
    if("RI" %in% type)
        res <- c(res, RI=sum(diag(nxx))/sum(nxx))

    if("J" %in% type)
        res <- c(res, J=nxx[2,2]/sum(nxx[-1]))
    
    if("FM" %in% type){
        tab <- table(x)
        w <- sum(tab*(tab-1))/2
        tab <- table(y)
        w <- w*sum(tab*(tab-1))/2
        res <- c(res, FM=nxx[2,2]/sqrt(w))
    }
    res        
})


countPairs <- function(x, y)
{
    if(length(x)!=length(y))
        stop("x an y must have the same length")
    
    res <- .C("countPairs",
              as.integer(x),
              as.integer(y),
              as.integer(length(x)),
              res=double(4))[["res"]]
    matrix(res, nrow=2, dimnames=list(0:1,0:1))
}
