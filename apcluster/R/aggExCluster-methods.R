aggExCluster.matrix <- function(s, x, includeSim=FALSE)
{
    noPriorClustering <- (missing(x) || is.null(x))

    if (length(dim(s)) != 2 || (ncol(s) != nrow(s) && noPriorClustering))
        stop("'s' must be a square matrix")

    AggResultObj <- new("AggExResult")

    K <- nrow(s)

    AggResultObj@l <- K

    preserveNames <- (length(rownames(s)) == nrow(s))

    if (noPriorClustering) ## no prior clustering
    {
        AggResultObj@maxNoClusters <- K
        AggResultObj@clusters[[K]] <- as.list(1:K)
        AggResultObj@exemplars[[K]] <- 1:K

        if (preserveNames)
        {
            AggResultObj@labels <- rownames(s)
            names(AggResultObj@exemplars[[K]]) <- rownames(s)

            for (i in 1:K)
                names(AggResultObj@clusters[[K]][[i]]) <- rownames(s)[i]
        }
        else
            AggResultObj@labels <- as.character(1:K)
    }
    else ## prior clustering
    {
        if (x@l != nrow(s))
            stop("data set sizes of 's' and 'x' do not match")

        AggResultObj@sel <- x@sel

        K <- length(x@exemplars)

        if (K < 1)
            stop("'x' empty or corrupted")

        AggResultObj@maxNoClusters <- K
        AggResultObj@clusters[[K]] <- x@clusters
        AggResultObj@exemplars[[K]] <- x@exemplars
        AggResultObj@labels <- paste("Cluster", 1:K)
    }

    if (K < 2)
    {
        warning("there is nothing to cluster")
        return(invisible(AggResultObj))
    }

    if (length(AggResultObj@sel) > 0)
    {
        colInd <- c(rep(0,nrow(s)))
        for (i in 1:length(AggResultObj@sel))
            colInd[AggResultObj@sel[i]] <- i
    }

    objMat <- matrix(NA, K, K) ## matrix of objective values for pairs
    exeMat <- matrix(NA, K, K) ## matrix of joint exemplars
    ## note: only the upper triangle of these matrices is non-NA

    actClust <- AggResultObj@clusters[[K]]
    actExem <- AggResultObj@exemplars[[K]]
    actLabels <- -(1:K)

    AggResultObj@merge <- matrix(NA, K - 1, 2)
    AggResultObj@height <- rep(0, K - 1)

    ## compute complete matrices before starting joining
    for (i in 1:(K - 1))
    {
        for (j in (i + 1):K)
        {
            joint <- c(actClust[[i]], actClust[[j]])

            if (length(AggResultObj@sel) > 0)
            {
                ci <- colInd[intersect(AggResultObj@sel,joint)]
                if (length(ci) > 0)
                {
                    cM <- colMeans(s[joint,
                                     colInd[intersect(AggResultObj@sel,joint)],
                                     drop=FALSE])
                    ex <- intersect(AggResultObj@sel,joint)[which.max(cM)]
                    exeMat[i, j] <- ex
                    objMat[i, j] <-
                        (mean(s[ex, colInd[intersect(AggResultObj@sel,
                                                     actClust[[i]])]]) +
                         mean(s[ex, colInd[intersect(AggResultObj@sel,
                                                     actClust[[j]])]])) / 2
                }
                else
                {
                    ## joining not possible - no similarities available
                    # exeMat[i,j] <- 0
                    # objMat[i,j] <- -Inf
                    stop("clusters cannot be joined because of missing ",
                         "similarity values;\n       maybe increasing the ",
                         "cluster size through decreasing\n",
                         "       the self similarity 'p' helps.")
                }
            }
            else
            {
                cM <- colMeans(s[joint, joint, drop=FALSE])
                ex <- joint[which.max(cM)]
                exeMat[i, j] <- ex
                objMat[i, j] <- (mean(s[ex, actClust[[i]]]) +
                                 mean(s[ex, actClust[[j]]])) / 2
            }

        }
    }

    ## agglomeration loop
    for (k in (K - 1):1)
    {
        tojoin <- which.max(objMat) - 1 ## determine pair to join
        I <- tojoin %% K + 1
        J <- floor(tojoin / K) + 1

        newClust <- c(actClust[[I]], actClust[[J]]) ## join them
        actClust[c(I, J)] <- NULL
        actClust[[k]] <- newClust
        actExem <- c(actExem[c(-I, -J)], exeMat[I, J])

        AggResultObj@clusters[[k]] <- actClust
        AggResultObj@merge[K - k, ] <- c(actLabels[I], actLabels[J])
        actLabels <- c(actLabels[c(-I, -J)], K - k)
        AggResultObj@height[K - k] <- objMat[I, J]
        AggResultObj@exemplars[[k]] <- actExem

        if (preserveNames)
            names(AggResultObj@exemplars[[k]]) <- colnames(s)[actExem]

        if (k == 1) break

        ## rearrange matrices objMat and exeMat
        ## put values for unchanged clusters in the first k-1 rows/columns
        indexVec <- 1:(k + 1)
        indexVec <- indexVec[c(-I, -J)]

        exeMat[1:(k - 1), 1:(k - 1)] <- exeMat[indexVec, indexVec, drop=FALSE]
        objMat[1:(k - 1), 1:(k - 1)] <- objMat[indexVec, indexVec, drop=FALSE]

        ## wipe out k+1-st column
        exeMat[, k + 1] <- NA
        objMat[, k + 1] <- NA

        ## update k-th column with objective values and joint exemplars of
        ## unchanged clusters and the newly joined cluster
        for (i in 1:(k - 1))
        {
            joint <- c(actClust[[i]], actClust[[k]])

            if (length(AggResultObj@sel) > 0)
            {
                ci <- colInd[intersect(AggResultObj@sel,joint)]
                if (length(ci) > 0)
                {
                    cM <- colMeans(s[joint,
                                     colInd[intersect(AggResultObj@sel,joint)],
                                     drop=FALSE])
                    ex <- intersect(AggResultObj@sel,joint)[which.max(cM)]
                    exeMat[i, k] <- ex
                    objMat[i, k] <- (mean(s[ex,
                                            colInd[intersect(AggResultObj@sel,
                                                             actClust[[i]])]]) +
                                     mean(s[ex,
                                            colInd[intersect(AggResultObj@sel,
                                                             actClust[[k]])]])) / 2
                }
                else
                {
                    ## joining not possible - no similarities available
                    # exeMat[i,j] <- 0
                    # objMat[i,j] <- -Inf
                    stop("clusters cannot be joined because of missing ",
                         "similarity values")
                }
            }
            else
            {
                cM <- colMeans(s[joint, joint, drop=FALSE])
                ex <- joint[which.max(cM)]
                exeMat[i, k] <- ex
                objMat[i, k] <- (mean(s[ex, actClust[[i]]]) +
                                 mean(s[ex, actClust[[k]]])) / 2
            }
        }
    }

    ## finally, determine reordering for dendrogram plotting
    AggResultObj@order <- determineOrder(AggResultObj@merge,
                                         AggResultObj@height, K - 1)

    AggResultObj@call <- deparse(sys.call(-1))

    if (includeSim)
        AggResultObj@sim <- s

    AggResultObj
}

setMethod("aggExCluster", signature("matrix", "missing" ), aggExCluster.matrix)
setMethod("aggExCluster", signature("matrix", "ExClust" ), aggExCluster.matrix)


aggExCluster.Matrix <- function(s, x, includeSim=FALSE)
{
    if (is(s, "sparseMatrix"))
    {
        s <- as.SparseSimilarityMatrix(s)

        rng <- range(s@x)

        fill <- 2 * rng[1] - rng[2]

        s <- as.DenseSimilarityMatrix(s, fill=fill)
    }
    else
        s <- as.DenseSimilarityMatrix(s)

    if (missing(x))
        res <- aggExCluster(s, includeSim=includeSim)
    else
        res <- aggExCluster(s, x, includeSim=includeSim)

    res
}

setMethod("aggExCluster", signature("Matrix", "missing" ), aggExCluster.Matrix)
setMethod("aggExCluster", signature("Matrix", "ExClust" ), aggExCluster.Matrix)


aggExCluster.Clust <- function(s, x, includeSim=TRUE)
{
    if (all(dim(x@sim) <= 1))
        stop("similarity matrix not included in object")

    AggResultObj <- aggExCluster(x@sim, x)

    AggResultObj@call <- deparse(sys.call(-1))

    if (includeSim)
        AggResultObj@sim <- x@sim

    AggResultObj
}

setMethod("aggExCluster", signature("missing" , "ExClust" ), aggExCluster.Clust)


aggExCluster.function <- function(s, x, includeSim=TRUE, ...)
{
    if (is.data.frame(x))
        x <- as.matrix(x[, sapply(x, is.numeric)])

    if (is.matrix(x))
        N <- nrow(x)
    else
        N <- length(x)

    if (N < 2) stop("cannot cluster less than 2 samples")

    if (!is.function(s))
    {
        if (!is.character(s) || !exists(s, mode="function"))
            stop("invalid distance function")

        s <- match.fun(s)
    }

    sim <- s(x=x, ...)

    if (!is.matrix(sim) || (nrow(sim) != N) || ncol(sim) != N)
        stop("computation of similarity matrix failed")

    AggResultObj <- aggExCluster(sim)

    AggResultObj@call <- deparse(sys.call(-1))

    if (includeSim)
        AggResultObj@sim <- sim

    AggResultObj
}

setMethod("aggExCluster", signature("function" , "ANY"), aggExCluster.function)
setMethod("aggExCluster", signature("character", "ANY"), aggExCluster.function)


## auxiliary function for determining the order for dendrogram plotting
## fills up order recursively starting from the last merge
determineOrder <- function(merge, height, k)
{
    I <- merge[k, 1] ## I and J are the clusters merged in the k-th step
    J <- merge[k, 2]

    if (I < 0 && J < 0) ## if both are singletons, list I first
        return(c(-I, -J))
    else if (I < 0) ## if I is a singleton and J is not, list it first
        return(c(-I, determineOrder(merge, height, J)))
    else if (J < 0) ## if J is a singleton and I is not, list it first
        return(c(-J, determineOrder(merge, height, I)))
    else ## if both are non-singleton clusters, list the "tighter" cluster
    {    ## on the left-hand side (see ?hclust)
        if (height[I] > height[J])
            return(c(determineOrder(merge, height, I),
                     determineOrder(merge, height, J)))
        else
            return(c(determineOrder(merge, height, J),
                     determineOrder(merge, height, I)))
    }
}
