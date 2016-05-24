##' Depends on the settings of adpclust, draw figures showing silhouette vs. number of clusters, f vs. delta with selected centroids, and original data (projected to the first two principal components if dim > 2) colored by cluster assignments. 
##'
##' @title Visualize the result of adpclust()
##' @param x an object of class "adpclust". Result of adpclust().
##' @param cols vector of colors used to distinguish different clusters. Recycled if necessary.
##' @param ... Not used.
##' @return NULL
##'
##' @export
##'
##' @examples
##' ## Load a data set with 3 clusters
##' data(clust3)
##' ## Automatically select cluster centroids
##' ans <- adpclust(clust3, centroids = "auto")
##' plot(ans)

plot.adpclust <- function(x,
                          cols = "default",
                          ...)
{
    if(cols == "default") cols = defCol()
    pars <- x$testPars
    ## if number of colors is not large enough, recycle
    nclust <- length(x$centers)
    if((temp <- ceiling(nclust / length(cols))) > 1)
        cols <- rep(cols, temp)[1:nclust]
    sil.vs.nclust <- score.vs.h <- FALSE
    f.vs.delta <- clust.res <- TRUE
    if(pars$centroids == "auto"){
        if(length(pars$nclusts) > 1)
            sil.vs.nclust <- TRUE
        else{
            if(length(pars$hs) > 1)
                score.vs.h <- TRUE
            else{}
        }
    }else{}

    par(mfrow = c(1, 2 + sil.vs.nclust + score.vs.h))

    ##--------------------
    ## nclust vs. silhouette
    ##--------------------    
    if(sil.vs.nclust){
        plot(pars$nclusts, pars$nclustScores, type = "b", xlab = "number of clusters",
             ylab = "silhouette", xaxt = "n", 
             main = "#clust vs silouette \n choosing nclust")
        abline(v = nclust, col = "red", lty = 2)
        axis(1, at = pars$nclusts, labels = pars$nclusts)
    }

    ##--------------------
    ## f vs delta (best one)
    ##--------------------
    if(f.vs.delta){ ## 
        plot(x$f, x$delta, xlab = "f(x)", ylab = "delta(x)", 
             main = "f(x) vs delta(x) \n chosen centers")
        f.range <- range(x$f)
        delta.range <- range(x$delta)
        points(x$f[x$centers], x$delta[x$centers], 
               col = cols, pch = 19, cex = 1.1)
        if(pars$centroids == "auto"){
            if(length(pars$f.cut) > 0){
                ## abline(v = f.range[1] + pars$f.cut * 
                ##            (f.range[2] - f.range[1]), col = "red", lty = 2)
                abline(v = quantile(x$f, probs = pars$f.cut), col = "red", lty = 2)
            }else{
                segments(f.range[1], delta.range[2],
                         f.range[2], delta.range[1], col = "red", lty = 2)
            }
        }
        text(x$f[x$center], x$delta[x$center], labels = x$center, cex = 0.6, pos = 1)

    }
    
    ##--------------------
    ## score vs. h
    ##--------------------
    if(score.vs.h){
        plot(score ~ h, pars$hs.scores, xlab = "h's", ylab = "silouette", 
             main = "selection with the best score")
        abline(v = x$h, col = "red", lty = 2)
    }
    
    ##--------------------
    ## clustering results. first 2 PC if dim > 2.
    ##--------------------
    maintitle <- ifelse(attr(x$dat, "type") == "pc", "Clustering Result \nFirst two principal components", "Clustering Result")
    myxlab <- ifelse(attr(x$dat, "type") == "pc", "PC 1", "x1")
    myylab <- ifelse(attr(x$dat, "type") == "pc", "PC 2", "x2")
    plot(x$dat, col = cols[x$cluster],
         main = maintitle,
         xlab = myxlab, ylab = myylab, pch = 1)
    points(x$dat[x$centers,], col = cols, pch = 4, cex = 3)
    ## legend("topright",
    ##        legend = c(paste0("cluster", 1:nclust)),
    ##        pch = 19,
    ##        col = c(cols[(1:nclust - 1) %% length(cols) + 1]),
    ##        cex = 0.6
    ##    )
    ## if(!is.null(file)) dev.off()
}

