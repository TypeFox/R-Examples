#
#  Copyright (C) 2009 Friedrich Leisch
#  $Id: family.R 3 2013-06-12 10:06:43Z leisch $
#

kccaFamily <- function(which=NULL, dist=NULL,
                       cent=NULL, name=which,
                       preproc=NULL,
                       trim=0, groupFun="minSumClusters")
{
    if(is.null(which) && is.null(dist))
        stop("either ", sQuote("which")," or ", sQuote("dist"),
             " must be specified\n")

    if(is.null(name)) name <- deparse(substitute(dist))

    z <- new("kccaFamily", name=name)

    if(!is.null(preproc)) z@preproc <- preproc

    if(!is.null(which)){
        which <- match.arg(which, c("kmeans", "kmedians",
                                    "angle", "jaccard",
                                    "ejaccard"))
        if(!is.null(name)) z@name <- which

        if(which == "kmeans"){
            z@dist <- distEuclidean
            if(trim==0){
                z@cent <- function(x) colMeans(x)
                z@wcent <- function(x, weights)
                    colMeans(x*normWeights(weights))
            }
            else{
                z@cent <- function(x)
                    apply(x, 2, mean, trim=trim)
                z@wcent <- function(x, weights)
                    apply(x*normWeights(weights), 2, mean, trim=trim)
            }
            z@weighted <- TRUE
        }
        else if(which == "kmedians"){
            z@dist <- distManhattan
            z@cent <- function(x) apply(x, 2, median)
        }
        else if(which == "angle"){
            z@dist <- distAngle
            z@cent <- centAngle
            z@wcent <- wcentAngle
            z@weighted <- TRUE
            z@preproc <- function(x) x/sqrt(rowSums(x^2))
        }
        else if(which == "jaccard"){
            z@dist <- distJaccard
            z@cent <- function(x) centOptim01(x, dist=distJaccard)
        }
        else if(which == "ejaccard"){
            z@dist <- distJaccard
            z@cent <- function(x) colMeans(x)
        }
    }
    else{
        if(is.character(dist))
            z@dist <- get(dist, mode="function")
        else
            z@dist <- dist

        if(is.null(cent))
            z@cent <- function(x)
                centOptim(x, dist=dist)
        else if(is.character(cent))
            z@cent <- get(cent, mode="function")
        else
            z@cent <- cent        
    }

    ## fill in z@cluster and z@allcent which use lexical scoping
    eval(FAMILY_CLUSTER_ALLCENT)
    
    if(is.character(groupFun))
        groupFun <- get(groupFun, mode="function")
    z@groupFun <- groupFun

    z
}


FAMILY_CLUSTER_ALLCENT <- expression({
    
    z@cluster <- function(x, centers, n=1, distmat=NULL){

        if(is.null(distmat))
            distmat <- z@dist(x, centers)
        
        if(n==1){
            return(max.col(-distmat))
        }
        else{
            r <- t(matrix(apply(distmat, 1,
                                rank, ties.method="random"),
                          nrow=ncol(distmat)))
            z <- list()
            for(k in 1:n)
                z[[k]] <- apply(r, 1, function(x) which(x==k))
        }
        return(z)
    }

    z@allcent <- function(x, cluster, k=max(cluster, na.rm=TRUE))
    {
        centers <- matrix(NA, nrow=k, ncol=ncol(x))
        for(n in 1:k){
            if(sum(cluster==n, na.rm=TRUE)>0){
                centers[n,] <- z@cent(x[cluster==n,,drop=FALSE])
            }
        }
        centers
    }
})

