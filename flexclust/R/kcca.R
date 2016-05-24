#
#  Copyright (C) 2005-2009 Friedrich Leisch
#  $Id: kcca.R 3 2013-06-12 10:06:43Z leisch $
#

normWeights <- function(x) x/mean(x)

mlogit <- function(x){
    x <- exp(-x)
    x/rowSums(x)
}
    
###**********************************************************


kcca <- function(x, k, family=kccaFamily("kmeans"), weights=NULL,
                 group=NULL, control=NULL, simple=FALSE, save.data=FALSE)
{
    MYCALL <- match.call()
    control <- as(control, "flexclustControl")
    x <- as(x, "matrix")
    x <- family@preproc(x)
    N <- nrow(x)
    
    if(control@classify=="auto"){
        control@classify="hard"
    }

    if(!is.null(group))
    {       
        if(length(group)>N)
            warning("group vector longer than nrow(x)")
        
        group <- rep(group, length=N)
        group <- as.factor(group)
    }
    
    if(!is.null(weights))
    {
        control@classify="weighted"

        if(!family@weighted)
            stop("Centroid function of family cannot use weights")

        if(!is.null(group))
            stop("Weighted clustering for grouped data is not implemented yet.")
        ## recycle to number of observations
        weights <- rep(weights, length=N)
    }
    
    centers <- initCenters(x, k, family, control)
    cluster <- centers$cluster
    k <- centers$k
    centers <- centers$centers

    sannprob <- control@simann[1]
    if(control@classify %in% c("hard", "simann"))
    {
        for(iter in 1:control@iter.max){
            
            clustold <- cluster
            distmat <- family@dist(x, centers)
            cluster <- family@cluster(x, distmat=distmat)

            if(control@classify == "simann"){
                cluster <- perturbClusters(cluster, sannprob)
                sannprob <- sannprob * control@simann[2]
            }

            if(!is.null(group))
                cluster <- family@groupFun(cluster, group, distmat)
            
            centers <- family@allcent(x, cluster=cluster, k=k)

            ## NAs in centers are empty clusters
            centers <- centers[complete.cases(centers),,drop=FALSE]
            k <- nrow(centers)
            
            changes <- sum(cluster!=clustold)
            if(control@verbose && (iter%%control@verbose==0)){
                td <- sum(distmat[cbind(1:N, cluster)])
                printIter(iter, paste(changes, format(td), sep=" / "),
                          "Changes / Distsum")
            }

            if(changes==0) break
        }
    }
    else if(control@classify=="weighted")
    {
        td <- -1
        for(iter in 1:control@iter.max)
        {
            td.old <- td
            distmat <- family@dist(x, centers)
            cluster <- family@cluster(distmat=distmat)

            td <- sum(distmat[cbind(1:N, cluster)])
            if(control@verbose && (iter%%control@verbose==0))
                printIter(iter, td, "Sum of distances")
                                
            if(abs(td-td.old)/(abs(td)+0.1) < control@tolerance) break

            ## for weight computation we need similarities
            distmat <- mlogit(distmat)

            for(n in 1:k){
                w <- weights*distmat[,n]
                w[cluster==n] <- w[cluster==n]+control@gamma
                centers[n,] <- family@wcent(x, weights=w)
            }
        }
        ## make sure the final centers are centroids of their partitions
        for(n in 1:k){
            centers[n,] <- family@cent(x[cluster==n,,drop=FALSE])
        }
    }
    else
        stop("Unknown classification method")

    centers <- centers[complete.cases(centers),,drop=FALSE]
    
    z <- newKccaObject(x=x, family=family, centers=centers, group=group,
                       iter=iter,
                       converged=(iter<control@iter.max),
                       call=MYCALL,
                       control=control,
                       simple=simple)

    if(save.data)
        z@data <- ModelEnvMatrix(designMatrix=x)

    z
}

###**********************************************************

initCenters <- function(x, k, family, control)
{
    cluster <- integer(nrow(x))
    if(is.matrix(k)){
        centers <- k
        if(any(duplicated(centers)))
            stop("initial centers are not distinct")
        k <- nrow(k)
    }
    else{
        if(length(k)>1){
            cluster <- rep(k, length=nrow(x))
            k <- max(cluster)
            centers <- family@allcent(x, cluster=cluster, k=k)
        }
        else{
            k <- as.integer(k)
            if(k<2) stop("number of clusters must be at least 2")
            ## we need to avoid duplicates here
            x <- na.omit(unique(x))
            if(nrow(x) < k)
                stop("k larger than number of distinct complete data points in x")
            centers <- do.call(control@initcent,
                               list(x=x, k=k, family=family))
        }
    }
    list(centers=centers, cluster=cluster, k=k)
}

randomcent <- function(x, k, family)
{
    x[sample(1:nrow(x), k), , drop=FALSE]
}

kmeanspp <- function(x, k, family)
{
    centers <- matrix(0, nrow=k, ncol=ncol(x))
    centers[1,] <- x[sample(1:nrow(x), 1), , drop=FALSE]
    d <- family@dist(x, centers[1L,,drop=FALSE])^2
    for(l in 2:k){
        centers[l,] <- x[sample(1:nrow(x), 1, prob=d), , drop=FALSE]
        d <- pmin(d, family@dist(x, centers[l,,drop=FALSE])^2)
    }
    centers
}
                  

###**********************************************************


newKccaObject <- function(x, family, centers, group=NULL, simple=FALSE,
                          ...)
{
    distmat <- family@dist(x, centers)

    z <- newKccasimpleObject(x=x, family=family, centers=centers, 
                             group=group, distmat=distmat, ...)
    
    if(!simple){
        z <- simple2kcca(x=x, from=z, group=group, distmat=distmat)
    }
    
    z
}

newKccasimpleObject <- function(x, family, centers, group=NULL,
                                distmat=NULL, ...)
{
    if(is.null(distmat))
        distmat <- family@dist(x, centers)

    cluster <- family@cluster(distmat=distmat)
    if(!is.null(group))
        cluster <- family@groupFun(cluster, group, distmat)
    names(cluster) <- rownames(x)
    colnames(centers) <- colnames(x)
    
    size <- as.vector(table(cluster))
    cldist <- as(distmat[cbind(1:nrow(x), cluster)], "matrix")
    cluster <- as(cluster, "integer")
    
    new("kccasimple",
        k=as(nrow(centers),"integer"),
        centers=centers,
        cluster=cluster,
        family=family,
        clusinfo=clusinfo(cluster, cldist, simple=TRUE),
        cldist=cldist,
        ...)
}

simple2kcca <- function(x, from, group=NULL, distmat=NULL)
{
    if(is.null(distmat))
        distmat <- from@family@dist(x, from@centers)

    cluster <- from@family@cluster(n=2, distmat=distmat)
    if(!is.null(group))
        cluster <- from@family@groupFun(cluster, group, distmat)
        
    if(ncol(distmat)>1){
        ## at least 2 clusters
        cldist <- cbind(distmat[cbind(1:nrow(x), cluster[[1]])],
                        distmat[cbind(1:nrow(x), cluster[[2]])])
        clsim <- computeClusterSim(distmat,cluster)
    }
    else{
        ## only one cluster
        cldist <- distmat
        clsim <- as.matrix(1)
    }
            
    xcent <- from@family@cent(x)
    totaldist <- sum(from@family@dist(x, matrix(xcent,nrow=1)))

    z <- as(from, "kcca")
    z@second <- as(cluster[[2]], "integer")
    z@xcent <- xcent
    z@xrange <- apply(x, 2, range, na.rm=TRUE)
    z@totaldist <- totaldist
    z@clsim <- clsim
    z@clusinfo <- clusinfo(cluster[[1]], cldist, simple=FALSE)
    z@cldist <- cldist

    z
}

###**********************************************************


clusinfo <- function(cluster, cldist, simple=FALSE)
### cluster: vector of cluster memberships
### cldist: matrix with 1 or 2 columns
{
    size <- as.vector(table(cluster))
    
    clusinfo <-
        data.frame(size=size,
                   av_dist = as.vector(tapply(cldist[,1], cluster, sum))/size)

    if(!simple){
        clusinfo <- cbind(clusinfo,
                          max_dist = as.vector(tapply(cldist[,1], cluster, max)),
                          separation = as.vector(tapply(cldist[,2], cluster, min)))
    }
    clusinfo
}

###**********************************************************

computeClusterSim <- function(distmat, cluster)
{
    K <- max(cluster[[1]])
    z <- matrix(0, ncol=K, nrow=K)

    for(k in 1:K){
        ok1 <- cluster[[1]]==k
        if(any(ok1)){
            for(n in 1:K){
                if(k!=n){
                    ok2 <- ok1 & cluster[[2]]==n
                    if(any(ok2)){
                        z[k,n] <- 2*sum(distmat[ok2,k]/
                                        (distmat[ok2,k]+distmat[ok2,n]))
                    }
                }
            }
            z[k,] <- z[k,]/sum(ok1)
        }
    }
    diag(z) <- 1
    z
}


###**********************************************************

# make random changes in cluster membership
perturbClusters <- function(cluster, prob)
{
    x <- runif(length(cluster))
    x <- (x < prob)
    if(any(x))
        cluster[x] <- sample(1:max(cluster), sum(x), replace=TRUE)

    cluster
}


###**********************************************************

setMethod("show", "kccasimple",
function(object)
{
    cat(class(object), "object of family", sQuote(object@family@name),"\n\n")
    cat("call:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    if(object@k<20){
        cat("\ncluster sizes:\n")
        print(table(clusters(object), useNA="ifany"))
    }
    else{
        cat("\n", sum(!is.na(object@cluster)), " points in ",
            object@k, " clusters", sep="")
        if(any(is.na(object@cluster))){
            cat(",", sum(is.na(object@cluster)), "outliers")
        }
        cat("\nDistribution of cluster sizes:\n", sep="")
        print(summary(as.integer(table(clusters(object)))))
    }
    cat("\n")
})
        
setMethod("summary", "kccasimple",
function(object)
{
    cat(class(object), "object of family", sQuote(object@family@name),"\n\n")
    cat("call:", deparse(object@call,0.75*getOption("width")),
        sep="\n")
    cat("\ncluster info:\n")
    print(object@clusinfo)
    cat("\n")
    if(!object@converged) cat("no ")
    cat("convergence after", object@iter, "iterations\n")
    if("distsum" %in% info(object))
        cat("sum of within cluster distances:", info(object, "distsum"),"\n")
})

###**********************************************************


setGeneric("neighbours",
function(object, symm=TRUE, ...) standardGeneric("neighbours"))

setMethod("neighbours", signature(object="kcca"),
function(object, symm=TRUE, ...)
{
    neighbours(list(object@cluster, object@second))
})

setMethod("neighbours", signature(object="list"),
function(object, symm=TRUE, ...)
{
    z <- list()
    for(n in 1:max(object[[1]], object[[2]])){
        if(symm)
            z[[n]] <- unique(c(object[[2]][object[[1]]==n],
                               object[[1]][object[[2]]==n]))
        else
            z[[n]] <- unique(object[[2]][object[[1]]==n])
    }
    z
})

###**********************************************************


setGeneric("ntable", function(object, ...)
           standardGeneric("ntable"))

setMethod("ntable", signature(object="kcca"),
function(object, freq=TRUE, distance=FALSE, ...)
{
    z <- table(object@cluster, object@second)

    if(distance){
        freq <- FALSE
        z[z==0] <- NA
        z <- length(object@cluster)-z
        diag(z) <- 0
    }
    
    if(!freq)
        z <- z/rowSums(z,na.rm=TRUE)

    z
})
        
    
###**********************************************************        

