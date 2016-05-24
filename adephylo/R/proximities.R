###########
# proxTips
###########
proxTips <- function(x, tips="all",
                      method=c("patristic","nNodes","oriAbouheif","Abouheif","sumDD"),
                     a=1, normalize=c("row","col","none"), symmetric=TRUE, useC=TRUE){

    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## handle arguments
    x <- as(x, "phylo4")
    method <- match.arg(method)
    normalize <- match.arg(normalize)
    N <- nTips(x)
    if(tips[1]=="all") { tips <- 1:N }
    tips <- getNode(x, tips)

    ## some checks
    if (is.character(checkval <- checkPhylo4(x))) stop(checkval)
    if(any(is.na(tips))) stop("wrong tips specified")

    ## compute distances
    distMethod <- method
    if(length(grep("Abouheif", distMethod)>1)){
        distMethod <- "Abouheif"
    }
    D <- distTips(x, tips=tips, method=distMethod, useC=useC)
    D <- as.matrix(D)

    ## compute proximities
    res <- (1/D)^a
    diag(res) <- 0

    ## handle Abouheif with diagonal (Abouheif1)
    if(method=="oriAbouheif"){
        sumMarg <- apply(res,1,sum)
        diag(res) <- (1-sumMarg)
        normalize <- "none" # not needed (already bistochastic)
        symmetric <- FALSE # not needed (aleady symmetric)
    }

    ## standardization
    if(normalize=="row") {
        res <- prop.table(res, 1)
    }

    if(normalize=="col") {
        res <- prop.table(res, 2)
    }

    ## re-symmetrize
    if(symmetric){
        res <- 0.5 * (res + t(res))
    }

    ## set the output
    return(res)

} # end proxTips
