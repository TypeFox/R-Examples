`ordihull.centroids` <- 
function (ordmod, ordiplot, groups, draw = c("lines", 
    "polygon"), show.groups, return.outliers = F, return.index = F, ...) 
{
    draw <- match.arg(draw)
	display <- "sites"
    pts <- scores(ordiplot, display = display, ...)
    if (!missing(show.groups)) {
        take <- groups %in% show.groups
        pts <- pts[take, , drop = FALSE]
        groups <- groups[take]
    }
    draw <- match.arg(draw)
    out <- seq(along = groups)
    inds <- names(table(groups))
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
    out <- seq(along = groups)
    inds <- names(table(groups))
    for (is in inds) {
        gr <- out[groups == is]
        if (length(gr) > 1) {
            X <- pts[gr, ]
            hpts <- grDevices::chull(X)
            hpts <- c(hpts, hpts[1])
            if (draw == "lines"){ 
                graphics::lines(X[hpts, ], ...)
            }else{ 
                graphics::polygon(X[hpts, ], ...)
            }
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

