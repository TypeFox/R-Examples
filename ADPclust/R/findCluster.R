#' Uses f and delta (calculated from findRd()) to find clusters.
#'
#' findCluster() plots f vs. delta plot, then allows user to select
#' centers of clusters and outliers. 
#' @title User-interactive Finding Clusters 
#' 
#' @param x a list. Return value from findRd2(). The list contains "f", "delta",
#' "dat", and "distm".
#' @param mycols a vector of character string. Colors to be used to distinguish
#' different cluster.
#' 
#' @return a list of the following items:
#' \itemize{
#' \item f: vector of f's.
#' \item delta: vector of delta's.
#' \item icenter: indices of cluster centers.
#' \item clusters: integer string recording cluster assignments. 
#' \item score: silouette of the clustering result.
#' \item h: bandwith h.
#' }
#' 
#' 
#' @examples
#' ## Internal examples

findCluster <- function(x, mycols = NULL)
{
    ##------------------------------------
    ## Color scheme; Simplify notations; Initiate res data frame
    ##------------------------------------
    if(is.null(mycols)) mycols <- defCol()

    f <- x$fd[[1]]$f
    delta <- x$fd[[1]]$delta
    dat <- x$dat
    distm <- x$distm
    
    ##------------------------------------
    ## Plot delta vs f. Click to select centers & outliers
    ##------------------------------------
    ## dev.new(width = 12, height = 6)
    ## par(mfrow = c(1,2), mar = c(7,4,4,3), mgp = c(3,1,0))
    par(mfrow = c(1,1), mar = c(7,4,4,3), mgp = c(3,1,0))
    plot(f, delta, xlab = "",
         ylab = paste0(expression(delta), "(x)"), main = "Decision Plot")
    mtext("f(x)\n Select centroids by left clicking \nPress 'ESC' to end selection", side = 1, line = 4)
    cat("Waiting user selection of centroids on the density-distance plot.\n")
    centers <- pickCenter(f, delta, col = mycols, labelcex = 0.6) # indices of centers
    frange <- range(f); drange <- range(delta)
    rect(xleft =  frange[1] + 0.3 * (frange[2] - frange[1]),
         ybottom = drange[1] + 0.4 * (drange[2] - drange[1]),
         xright = frange[1] + 0.7 * (frange[2] - frange[1]),
         ytop = drange[1] + 0.6 * (drange[2] - drange[1]), col = "white")
    text(mean(frange), mean(drange), labels = "Selection Finished.")
         
    if(length(centers) < 2) stop("Select at least two centers.")
    ## ioutliers <- pickOutlier(f, delta, col = "black") # indices of centers
    ## Mark selected points
    ## noutlier <- length(ioutliers)
    ## res$outlier <- rep(FALSE, length(f))
    ## res$outlier[ioutliers] <- TRUE

    ##------------------------------------
    ## Assign pts to clusters
    ##------------------------------------
    ## 'rdist' finds distances between pts to centers; "euclidean distance";
    dist.to.c <- x$distm[, centers]
    cluster <- apply(dist.to.c, 1, FUN = which.min)
    score <- mean(silhouette(cluster, x$distm)[,3])
    
    ## dist.to.c <- rdist(dat, dat[centers, ]) # length(dat) by length(centers) distance matrix.
    ## score <- intCriteria(traj = as.matrix(dat), part = res$cluster, crit = "Silhouette")[[1]]

    ##------------------------------------
    ## Find borders and cores (!halos)
    ##------------------------------------
    ## Find borders.
    ## distm.diffc <- distm / outer(res$cluster, res$cluster, FUN = "!=")
    ## res$on.border <- apply(distm.diffc, 2, min, na.rm = TRUE) < x$h
    ## ## Find cluster cores
    ## f1 <- f 
    ## f1[!res$on.border] <- -Inf # replacing !on.border values in f by -Inf
    ## res$core <- f > ave(f1, res$cluster, FUN = max)

    ##------------------------------------
    ## Return answers
    ##------------------------------------
    ans <- vector("list", 5)
    ans$cluster <- cluster
    ans$centers <- centers
    ans$score <- score
    ans$nclust <- length(centers)
    ans$h <- x$h
    class(ans) <- c("ADPclust", "list")
    invisible(ans)
}



