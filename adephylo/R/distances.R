###########
# distTips
###########
distTips <- function(x, tips="all",
                      method=c("patristic","nNodes","Abouheif","sumDD"), useC=TRUE){

    ## if(!require(phylobase)) stop("phylobase package is not installed")

    if(useC){
        tre <- as(x, "phylo")
        n <- as.integer(nTips(tre))
        resSize <- as.integer(n*(n-1)/2)
        res <- double(resSize)
        method <- match.arg(method)
        method <- match(method, c("patristic","nNodes","Abouheif","sumDD"))
        if(is.null(tre$edge.length)){
            tre$edge.length <- as.double(rep(1, nrow(tre$edge)))
        }

        temp <- .C("distalltips", as.integer(tre$edge[,1]), as.integer(tre$edge[,2]), as.double(tre$edge.length), nrow(tre$edge), n, res, resSize, as.integer(method), PACKAGE="adephylo")
        res <- temp[[6]]

        class(res) <- "dist"
        attr(res, "Size") <- nTips(tre)
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        attr(res, "method") <- paste("Phylogenetic: ",method,sep="")
        attr(res, "call") <- match.call()
        attr(res, "Labels") <- tre$tip.label
    } else {

        ## handle arguments
        x <- as(x, "phylo4")
        method <- match.arg(method)
        N <- nTips(x)
        if(tips[1]=="all") { tips <- 1:N }
        tips <- getNode(x, tips)
        tips.names <- names(tips)

        ## some checks
        if (is.character(checkval <- checkPhylo4(x))) stop(checkval)
        if(any(is.na(tips))) stop("wrong tips specified")

        ## create all couples of observations
        findAllPairs <- function(vec){
            res <- list(i=NULL,j=NULL)
            k <- 0
            for(i in 1:(length(vec)-1)){
                for(j in (i+1):length(vec)){
                    k <- k+1
                    res[[1]][k] <- i
                    res[[2]][k] <- j
                }
            }
            res <- data.frame(res)
            return(res)
        }

        allPairs <- findAllPairs(tips) # this contains all possible pairs of tips

        ## get the shortest path between all pairs of tips
        if(method != "patristic") {
            allPath <- sp.tips(x, allPairs$i, allPairs$j, useTipNames=TRUE, quiet=TRUE)
        } else {
            allPath <- sp.tips(x, allPairs$i, allPairs$j, useTipNames=TRUE, quiet=TRUE,
                               include.mrca=FALSE)
        }

        ## compute distances
        if(method=="patristic"){
            if(!hasEdgeLength(x)) stop("x does not have branch length")
            ## add tip1 and tip2 to the paths, so that these edges are counted
            allPath.names <- names(allPath)
            allPath <- lapply(1:length(allPath), function(i)
                              c(allPath[[i]], allPairs[i,1], allPairs[i,2]) )
            names(allPath) <- allPath.names

            edge.idx <- lapply(allPath, function(e) getEdge(x, e) ) # list of indices of edges
            allEdgeLength <- edgeLength(x)
            res <- lapply(edge.idx, function(idx) sum(allEdgeLength[idx], na.rm=TRUE) )
        } # end patristic

        if(method=="nNodes"){
            res <- lapply(allPath, length)
        } # end nNodes

        if(method=="Abouheif"){
            E <- x@edge
            f1 <- function(onePath){ # computes product of dd for one path
                temp <- table(E[,1])[as.character(onePath)] # number of dd per node
                return(prod(temp))
            }
            res <- lapply(allPath, f1)
        } # end Abouheif

        if(method=="sumDD"){
            E <- x@edge
            f1 <- function(onePath){ # computes sum of dd for one path
                temp <- table(E[,1])[as.character(onePath)] # number of dd per node
                return(sum(temp))
            }
            res <- lapply(allPath, f1)
        } # end sumDD

        ## convert res to a dist object
        res <- unlist(res)
        class(res) <- "dist"
        attr(res, "Size") <- length(tips)
        attr(res, "Diag") <- FALSE
        attr(res, "Upper") <- FALSE
        attr(res, "method") <- paste("Phylogenetic: ",method,sep="")
        attr(res, "call") <- match.call()
        attr(res, "Labels") <- tips.names
    }
    return(res)

} # end distTips







###########
# distRoot
###########
distRoot <- function(x, tips="all", method=c("patristic","nNodes","Abouheif","sumDD") ){
    ## if(!require(phylobase)) stop("phylobase package is not installed")

    ## handle arguments
    x <- as(x, "phylo4")
    method <- match.arg(method)
    N <- nTips(x)
    if(tips[1]=="all") { tips <- 1:N }
    tips <- getNode(x, tips)
    tips.names <- names(tips)
    x <- as(x, "phylo4")
    root <- getNode(x, N+1) # so that we have a named node

    ## some checks
    if(is.character(checkval <- checkPhylo4(x))) stop(checkval)
    if(any(is.na(tips))) stop("wrong tips specified")


    ## main computations

    ## get path from root to tops
    allPath <- lapply(tips, function(tip) .tipToRoot(x, tip, root, include.root = TRUE))

    ## compute distances
    if(method=="patristic"){
        if(!hasEdgeLength(x)) stop("x does not have branch length")
        ## add the concerned tips to the paths, so that these edges are counted
        allPath.names <- names(allPath)
        allPath <- lapply(1:length(allPath), function(i) c(allPath[[i]], tips[i]) )
        names(allPath) <- allPath.names

        edge.idx <- lapply(allPath, function(e) getEdge(x, e) ) # list of indices of edges
        allEdgeLength <- edgeLength(x)
        res <- sapply(edge.idx, function(idx) sum(allEdgeLength[idx], na.rm=TRUE) )
    } # end patristic

    if(method=="nNodes"){
        res <- sapply(allPath, length)
    } # end nNodes

    if(method=="Abouheif"){
        E <- x@edge
        f1 <- function(onePath){ # computes product of dd for one path
            temp <- table(E[,1])[as.character(onePath)] # number of dd per node
            return(prod(temp))
        }

        res <- sapply(allPath, f1)
    } # end Abouheif

    if(method=="sumDD"){
        E <- x@edge
        f1 <- function(onePath){ # computes sum of dd for one path
            temp <- table(E[,1])[as.character(onePath)] # number of dd per node
            return(sum(temp))
        }

        res <- sapply(allPath, f1)
    } # end sumDD


    ## the output is a named numeric vector
    return(res)
} # end distRoot
