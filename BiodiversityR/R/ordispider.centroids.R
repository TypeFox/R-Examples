`ordispider.centroids` <-
function (ordmod, ordiplot, groups, w = weights(ordiplot, 
    "sites"), show.groups, return.outliers = F, return.index = F, ...) 
{
    if (inherits(ordiplot, "cca") && missing(groups)) {
        lc <- scores(ordiplot, display = "lc", ...)
        wa <- scores(ordiplot, display = "wa", ...)
        graphics::segments(lc[, 1], lc[, 2], wa[, 1], wa[, 2], ...)
        return(invisible())
    }
	display <- "sites"
    ordispidercentres <- function(ordiplot, groups, display = "sites", 
        w = weights(ordiplot, display), show.groups) {
        pts <- scores(ordiplot, display = display)
        w <- eval(weights(ordiplot, display))
        if (length(w) == 1) 
            w <- rep(1, nrow(pts))
        if (is.null(w)) 
            w <- rep(1, nrow(pts))
        if (!missing(show.groups)) {
            take <- groups %in% show.groups
            pts <- pts[take, , drop = FALSE]
            groups <- groups[take]
            w <- w[take]
        }
        out <- seq(along = groups)
        inds <- names(table(groups))
        result <- array(NA, dim = c(length(inds), 2))
        rownames(result) <- inds
        for (is in inds) {
            gr <- out[groups == is]
            if (length(gr) > 1) {
                X <- pts[gr, ]
                W <- w[gr]
                ave <- apply(X, 2, weighted.mean, w = W)
                result[is, 1] <- ave[1]
                result[is, 2] <- ave[2]
            }
        }
        return(result)
    }
    pts <- scores(ordiplot, display = display, ...)
    w <- eval(w)
    if (length(w) == 1) 
        w <- rep(1, nrow(pts))
    if (is.null(w)) 
        w <- rep(1, nrow(pts))
    if (!missing(show.groups)) {
        take <- groups %in% show.groups
        pts <- pts[take, , drop = FALSE]
        groups <- groups[take]
        w <- w[take]
    }
    out <- seq(along = groups)
    inds <- names(table(groups))
    origcentres <- ordispidercentres(ordiplot = ordiplot, groups = groups)
    if (is.null(ordmod$CCA)) 
        stop(paste("function requires ordination result with constrained axes"))
    mod2 <- summary(ordmod)
    caxes <- length(mod2$ev.con)
    mod2 <- summary(ordmod,axes=caxes)
    scoord <- as.matrix(mod2$sites)[,1:caxes]
    out <- seq(along = groups)
    inds <- names(table(groups))
    g <- length(inds)
    centres <- array(dim=c(g,caxes),NA)
    rownames(centres) <- inds
    for (i in 1:g) {
        ind1 <- groups == inds[i]
        centroids <- apply(scoord[ind1,],2,mean)
        centroids <- t(as.matrix(centroids))
        for (j in 1:caxes) {
            centres[i,j] <- centroids[1,j] 
        }
    }
    index <- rep(FALSE, length(out))
    membership <- array(dim=c(length(out),4))
    colnames(membership) <- c("own.centroid","closest.centroid","own.centroid.dist","closest.centroid.dist")
    for (i in 1:length(index)) {
        mat1 <- rbind(centres, scoord[i, ])
        distances <- vegdist(mat1, "euc")
        distances <- as.matrix(distances)
        distances <- distances[1:g, g + 1]
        ind2 <- groups[i] == names(distances)
        if (distances[ind2] == min(distances)) {
            index[i] <- TRUE
        }
        membership[i,1] <- names(distances)[groups[i]]
        membership[i,2] <- names(distances)[which.min(distances)]
        membership[i,3] <- distances[groups[i]]
        membership[i,4] <- distances[which.min(distances)]
    }
    groups <- groups[index]
    pts <- pts[index, ]
    w <- w[index]
    out <- seq(along = groups)
    inds <- names(table(groups))
    for (is in inds) {
        gr <- out[groups == is]
        origrow <- rownames(origcentres)==is
        if (length(gr) > 1) {
            X <- pts[gr, ]
            W <- w[gr]
            ave1 <- origcentres[origrow,1]
            ave2 <- origcentres[origrow,2]
            graphics::segments(ave1, ave2, X[, 1], X[, 2], ...)
        }
    }
    if (return.outliers == T) {
        rownames(membership) <- rownames(scoord)
        index <- index==F
        membership <- membership[index,]
        membership <- data.frame(membership)
        return(membership)
    }
    if (return.index == T) {
        cat("Sites that are closest to their own centroid\n")
        return(index)
    }
    invisible()
}



