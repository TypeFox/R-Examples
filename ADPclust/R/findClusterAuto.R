#' @title Automaticly Find Clusters
#' 
#' @description This function uses rho and delta (calculated from findRd()) to
#' automaticly find clusters.
#' 
#' @details centers of clusters are selected as follows. Points corresponding to the largest nclust delta values are selected from
#' the points whose rho's are greater than 10th percentile. 
#' @param x a list. Return value from findRd(). The list contains "f", "delta",
#' "dat", and "distm".
#' @param nclust a number. Number of clusters.
#' @param ac an integer indicating which automatic cut method is used. Currently take two values:
#' \itemize{
#' \item{ac = 1: }{in the f vs. delta decision plot, 'nclust' points with f > percentile f.cut and nclust largest delta's are declaired centroids.}
#' \item{ac = 2: }{in the f vs. delta decision plot, denote by l the diagonal line connecting the point with smallest f and largest delta, and the point with largest f and smallest delta. 'nclust' points that are above l, and have are farthest away from l are declared centroids.}
#' }
#' @param mycols a vector of character string. Colors to be used to distinguish.
#' different cluster.
#' @param verbose if TRUE progress will be displayed.
#' @param findSil if FALSE silhouette score is NOT calculated, and the field that stores silhouette is set to -Inf. The default is TRUE. This argument is IGNORED if length(h) == 0 or length(h) > 1 or length(nclust) > 1, as silhouette scores are needed to select the best clustering in these cases.#' 
#' @param f.cut numeric vector. f.cut is used in variation = "auto" to automatically select cluster centroids from the decision plot. Points with f(x) > f.cut and high delta(x) are selected as candidate centroids. 
#' #' 
#' @return a list of the following items:
#' \itemize{
#' \item rho: vector of rho's.
#' \item delta: vector of delta's.
#' \item icenter: indices of cluster centers.
#' \item clusters: integer string recording cluster assignments. 
#' \item score: silouette of the clustering result.
#' \item h: bandwith h.
#' }
#' 
#' @examples
#' ## To be finished
findClusterAuto <- function(x,
                            mycols = NULL,
                            nclust = 2,
                            ac = 1,
                            f.cut = c(0.1, 0.2, 0.3),
                            verbose = FALSE,
                            findSil = TRUE)

{
    ##------------------------------------
    ## Color scheme
    ##------------------------------------
    if(is.null(mycols)) mycols <- defCol()
    h <- x$h; dat <- x$dat
    
    ##------------------------------------
    ## Find centers
    ##------------------------------------
    if(ac == 1){
        centers <- array(list(NULL), c(length(h), length(nclust), length(f.cut)),
                         dimnames = list(h = round(h,2), nclust = nclust, fCut = f.cut))
        for(i in 1:length(h)){
            f <- x$fd[[i]]$f
            delta <- x$fd[[i]]$delta
            for(k in 1:length(f.cut)){
                ##f0 <- min(f) + f.cut[j] * (max(f) - min(f))
                f0 <- quantile(f, probs = f.cut[k])
                delta1 <- delta
                delta1[f < f0] <- -Inf
                cts <- order(delta1, decreasing = TRUE)[1:max(nclust)]
                for(j in 1:length(nclust)){
                    centers[[i, j, k]] <- cts[1:nclust[j]]
                }
            }
        }
    }else if(ac == 2){
        centers <- array(list(NULL), c(length(h), length(nclust)),
                         dimnames = list(h = round(h, 2), nclust = nclust))
        for(i in 1:length(h)){
            f <- x$fd[[i]]$f
            delta <- x$fd[[i]]$delta
            x1 <- min(f); y1 <- max(delta)
            x2 <- max(f); y2 <- min(delta)
            pl.dist <- ((x2 - x1) * delta - (y2 - y1) * f + y2 * x1 - x2 * y1) /
                sqrt((y2 - y1) ^ 2 + (x2 - x1) ^ 2)
            cts <- order(pl.dist, decreasing = TRUE)[1:max(nclust)]
            for(j in 1:length(nclust)){
                centers[[i, j]] <- cts[1:nclust[j]]
            }
        }
    }else{
        stop("Wrong ac. Try ac = 1 or ac = 2.")
    }
    
    ##------------------------------------
    ## Find best scores from one cluster
    ##------------------------------------
    findScoreOneNclust <- function(clustID){
        if(verbose) cat("Working on #cluster =", nclust[clustID], "...")
        ## Extract Grid of centers of the current nclust
        if(ac == 1){
            myGrid <- centers[, clustID, ]
            hGrid <- cbind(1:length(h)) %*% rbind(rep(1, length(f.cut)))
            fCutGrid <- cbind(rep(1, length(h))) %*% rbind(1:length(f.cut))
            ##scores <- matrix(NA, nrow = nrow(myGrid), ncol = ncol(myGrid))
        }else{
            myGrid <- centers[, clustID]
            hGrid <- cbind(1:length(h)) %*% rbind(rep(1, length(f.cut)))
            fCutGrid <- NULL
            ##scores <- vector("numeric", length(myGrid))
        }
        ## Find unique centers
        nUniq <- 1
        unqCenters <- list()
        unqCenters[[1]] <- myGrid[[1]]
        attr(unqCenters[[1]], "h.id") <- 1
        attr(unqCenters[[1]], "f.cut.id") <- fCutGrid[[1]]
        if(length(myGrid) > 1){
            for(i in 2:length(myGrid)){
                loc <- belongsTo(myGrid[[i]], unqCenters)
                if(loc){ ## set i is NOT unique
                }else{ ## set i is unique
                    nUniq <- nUniq + 1
                    unqCenters[[nUniq]] <- myGrid[[i]]
                    attr(unqCenters[[nUniq]], "h.id") <- hGrid[[i]]
                    attr(unqCenters[[nUniq]], "f.cut.id") <- fCutGrid[[i]]
                }
            }
        }
        ## Find scores
        findOneScore <- function(oneSetCenters){
            dist.to.c <- x$distm[, oneSetCenters]
            clusters <- apply(dist.to.c, 1, FUN = which.min)
            if(findSil) score <- mean(silhouette(clusters, x$distm)[,3])
            else score = -Inf
            return(list(score = score, clusters = clusters))
        }
        ress <- lapply(unqCenters, findOneScore)
        scores <- sapply(ress, function(x) x$score)
        best <- which.max(scores)
        best.score <- scores[best]
        attr(best.score, "h") <- h[attr(unqCenters[[best]], "h.id")]
        attr(best.score, "h.id") <- attr(unqCenters[[best]], "h.id")
        attr(best.score, "f.cut") <- f.cut[attr(unqCenters[[best]], "f.cut.id")]
        attr(best.score, "nclust") <- nclust[[clustID]]
        h.score <- data.frame(h = h[sapply(unqCenters, function(x) attr(x, "h"))], scores = scores)
        hs.scores <- group_by(h.score, h) %>% dplyr::summarize(score = max(scores))
        if(verbose) cat("Done.\n")
        return(list(score = best.score,
                    cluster = ress[[best]]$clusters,
                    centers = unqCenters[[best]], 
                    hs.scores = hs.scores))
    }

    ##------------------------------------
    ## Main scripts
    ##------------------------------------
    scores.list <- lapply(1:length(nclust), FUN = findScoreOneNclust) # Best scores from each nclust
    scores <- sapply(scores.list, function(x) x$score) # extract just the scores
    best <- which.max(scores)
    
    ans <- list()
    ans$cluster <- scores.list[[best]]$cluster
    ans$centers <- scores.list[[best]]$centers
    ans$score <- scores[[best]]
    ans$nclust <- nclust[best]
    ans$h <- attr(scores.list[[best]]$score, "h")
    ans$f.cut <- attr(scores.list[[best]]$score, "f.cut")
    ans$hs.scores <- scores.list[[best]]$hs.scores
    ans$nclustScores <- scores
    ans$h.id <- attr(scores.list[[best]]$score, "h.id")
    invisible(ans)
}



