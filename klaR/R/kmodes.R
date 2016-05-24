kmodes <- function(data, modes, iter.max = 10, weighted = FALSE) {
### updates the mode of cluster "num",
### is called everytime an object is switches clusters
    update_mode <- function(num, num_var, data, cluster) {
        ## gather all objects of cluster "num"
        clust <- data[which(cluster == num),] 
        apply(clust, 2, function(cat) {
            ## compute the most frequent category for each variable
            cat <- table(cat)
            names(cat)[which.max(cat)]
        })
    } # end update mode


### computes the weighted distance between an object "obj" and the mode of its cluster
### the frequencies of categories in "data" are used for the weighting
    distance <- function(mode, obj, weights){
        if(is.null(weights))
            return(sum(mode != obj))
            
        obj <- as.integer(obj)
        different <- which(mode != obj) 
        n_mode <- n_obj <- numeric(length(different))
        for (i in seq(along = different)) { 
            ## frequencies are computed only if neccessary (not if distance is 0 anyway)
            weight <- weights[[different[i]]]
            names <- names(weight)
            n_mode[i] <- weight[which(names==mode[different[i]])]
            n_obj[i] <- weight[which(names==obj[different[i]])]
        }
        dist <- sum((n_mode + n_obj) / (n_mode * n_obj))
        return(dist)
    } # end distance


    n <- nrow(data)
    num_var <- ncol(data)
    data <- as.data.frame(data)
    
    cluster <- numeric(n) # the actual partition
    names(cluster) <- 1:n
    
    if (missing(modes)) 
        stop("'modes' must be a number or a matrix.")
    if (iter.max < 1) 
        stop("'iter.max' must be positive.")

    if (length(modes) == 1) { 
        ## if the desired number "k" of clusters is given, choose k object randomly from "data"
        k <- modes
        modes <- unique(data)[sample(nrow(unique(data)))[1:k],] 
        for (i in 1:k)
            cluster[which(rownames(data) == rownames(modes)[i])] <- i
    } else { 
        ## use given modes if appropriate
        if (any(duplicated(modes))) 
            stop("Initial modes are not distinct.")
        if (ncol(data) != ncol(modes)) 
            stop("'data' and 'modes' must have same number of columns")
        modes <- as.matrix(modes)
        k <- nrow(modes)
    }
    if (k > nrow(unique(data))) 
        stop("More cluster modes than distinct data points.")

    if(weighted){
        ## compute the frequencies of each category for each variable
        weights <- vector("list", num_var)
        for (i in 1:num_var) 
            weights[[i]] <- table(data[,i])
    } else {
        weights <- NULL
    }
    
    for (j in which(cluster==0)) { 
        ## first put all not yet assigned objects into the cluster, which has the nearest mode
        dist <- apply(modes, 1, function(x) distance(x, data[j,], weights))
        cluster[j] <- which.min(dist)
        modes[cluster[j],] <- update_mode(cluster[j], num_var, data, cluster) 
        ## update the mode of the cluster the object has been assigned to
    }
    for (i in 1:iter.max) {
        continue <- FALSE
        for (j in 1:n) { 
            ## run through all objects and assign them to the cluster with the nearest mode
            ## (or leave everything as is, if the object's current cluster's mode is still the nearest)
            dist <- apply(modes, 1, function(x) distance(x, data[j,], weights))
            clust_new <- which.min(dist)
            clust_old <- cluster[j]
            if (clust_new != clust_old) { 
                ## update the modes of old and new cluster, if object has switched from "clust_old" to "clust_new"
                    cluster[j] <- clust_new
                modes[clust_new,] <- update_mode(clust_new, num_var, data, cluster)
                modes[clust_old,] <- update_mode(clust_old, num_var, data, cluster)
                continue <- TRUE
            }
        }
        ## break, if during a complete iteration no object has changed clusters
        if (!continue) break
    } # end iterations

    cluster.size <- table(cluster) # sizes of all clusters
    if (length(cluster.size) < k) 
        warning("One or more clusters are empty.")
    ## compute for each cluster the sum of the distances of all objects of a cluster to the cluster's mode:
    diffs <- numeric(k)
    for (i in seq_along(cluster.size)) 
        diffs[i] <- sum(apply(data[cluster == i,], 1, function(x) sum(x != modes[i,]) ))
    rownames(modes) <- 1:k
    colnames(modes) <- colnames(data)
    result <- list(cluster = cluster, size = cluster.size, modes = modes, 
        withindiff = diffs, iterations = i, weighted = weighted)
    class(result) <- "kmodes"
    return(result)
}


print.kmodes <- function(x, ...)
{
    cat("K-modes clustering with ", length(x$size), " clusters of sizes ",
        paste(x$size, collapse=", "), "\n", sep="")
    cat("\nCluster modes:\n")
    print(x$modes, ...)
    cat("\nClustering vector:\n")
    print(x$cluster, ...)
    cat("\nWithin cluster simple-matching distance by cluster:\n")
    print(x$withindiff, ...)
    cat("\nAvailable components:\n")
    print(names(x))
    invisible(x)
}
