#
#  Copyright (C) 2006 Friedrich Leisch
#  $Id: grouping.R 3 2013-06-12 10:06:43Z leisch $
#


## Assign each observation to the cluster minimizing the sum
## of distances to all group members.
minSumClusters <- function(cluster, group, distmat)
{
    G <- levels(group)
    x <- matrix(0, ncol=ncol(distmat), nrow=length(G))

    for(n in 1:length(G)){
        x[n,] <- colSums(distmat[group==G[n],,drop=FALSE])
    }

    m <- max.col(-x)
    names(m) <- G
    z <- m[group]
    names(z) <- NULL

    if(is.list(cluster))
    {
        ## get second best
        x[cbind(1:nrow(x), m)] <- Inf
        m <- max.col(-x)
        names(m) <- G
        z1 <- m[group]
        names(z1) <- NULL
        z <- list(z, z1)
    }
    z
}
    


## Assign each observation to the cluster where the majority from its
## group belongs to
majorityClusters <- function(cluster, group, distmat)
{
    if(is.list(cluster))
    {
        K <- max(unlist(cluster))
        ## use factors to make sure all possible clusters appear in
        ## both tables
        x <- 2*table(group, factor(cluster[[1]], levels=1:K)) +
               table(group, factor(cluster[[2]], levels=1:K))
    }
    else{
        x <- table(group, cluster)
    }

    m <- max.col(x)
    names(m) <- row.names(x)
    z <- m[group]
    names(z) <- NULL
    
    if(is.list(cluster))
    {
        ## get second best
        x[cbind(1:nrow(x), m)] <- 0
        m <- max.col(x)
        names(m) <- row.names(x)
        z1 <- m[group]
        names(z1) <- NULL
        z <- list(z, z1)
    }
    z
}

    
## Assign each observation to a cluster where no other member from its
## group belongs to
differentClusters <- function(cluster, group, distmat)
{
    if(max(table(group)) > ncol(distmat))
        stop("Number of group members must be smaller or equal as number of clusters")
       
    INF <- 2*sum(distmat, na.rm=TRUE)
    distmat[is.na(distmat)] <- INF
    if(is.list(cluster)){
        z <- getDifferentCluster(cluster[[1]], group, distmat)
        distmat[cbind(1:nrow(distmat), z)] <- INF
        cluster <- list(z,
                        getDifferentCluster(cluster[[2]], group, distmat))
    }
    else{
        cluster <- getDifferentCluster(cluster, group, distmat)
    }

    cluster    
}
    
getDifferentCluster <- function(cluster, group, distmat)
{
    x <- table(group, cluster)
    ok <- (apply(x, 1, max)==1)
    nok.names <- unique(row.names(x[!ok,,drop=FALSE]))
    require("clue")

    for(n in nok.names){
        ok <- group==n
        if(sum(ok)>1)
            cluster[ok] <- solve_LSAP(distmat[ok,])
    }
    cluster
}

## solve_LSAP1 <- function (x, maximum = FALSE) 
## {
##     require("clue")
##     if (!is.matrix(x) || any(x < 0)) 
##         stop("x must be a matrix with nonnegative entries.")

##     if(nrow(x)>ncol(x))
##         stop("x must have less or equal rows than columns")

##     nr <- nrow(x)
##     nc <- ncol(x)
##     if(ncol(x) > nrow(x))
##         x <- rbind(x, matrix(2*sum(x), nrow=(ncol(x)-nrow(x)), ncol=ncol(x)))

##     if (maximum)
##         x <- max(x) - x
    
##     storage.mode(x) <- "double"
##     out <- .C("solve_LSAP", x, nc, p = integer(nc), PACKAGE = "clue")$p + 1
##     out <- out[1:nr]
##     class(out) <- "solve_LSAP"
##     out
## }    
