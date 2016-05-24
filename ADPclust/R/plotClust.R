plotClust <- function(dat, cres, file, mycols, plotType,
                      nclust, scores, criterion, centroids, ...)
{
    if(is.null(mycols)) mycols <- defCol()
    
    if(!is.null(file)){
        defaults <- list(filename = file, res = 300,
                         units = "in", width = 12, height = 4)
        args <- modifyList(defaults, list(...))
        do.call("png", args)
        ## png(filename = file, res = 300,
        ##     units = "in", width = 12, height = 4, ...)
    }

    ## If number of colors is not high enough, recycle
    Final.nclust <- length(cres$icenters)
    if((temp <- ceiling(Final.nclust / length(mycols))) > 1)
        mycols <- rep(mycols, temp)[1:Final.nclust]

    if(plotType == "twoplots"){
        par(mfrow = c(1,2))
    }else if(plotType %in% c("threeplots", "threeplots2")){
        par(mfrow = c(1,3))
    }

    ##--------------------
    ## nclust vs. silhouette
    ##--------------------    
    if(plotType == "threeplots"){
        plot(nclust, scores, type = "b", xlab = "# Cluster",
             ylab = "Silhouette", xaxt = "n", 
             main = "#clust vs silouette \n Choosing nclust")
        abline(v = Final.nclust, col = "red", lty = 2)
        axis(1, at = nclust, labels = nclust)
    }

    ##--------------------
    ## Final f vs delta
    ##--------------------
    plot(cres$f, cres$delta, xlab = "f(x)", ylab = "delta(x)", 
         main = "f(x) vs delta(x) \n chosen centers")
    f.range <- range(cres$f)
    delta.range <- range(cres$delta)
    if(!is.null(cres$cut.type)){
        if(cres$cut.type == 1){
            abline(v = f.range[1] + cres$cutvalue * 
                       (f.range[2] - f.range[1]), col = "red", lty = 2)
        }else{
            ## minf.loc <- which.min(cres$f)
            ## x1 <- cres$f[minf.loc]; y1 <- cres$delta[minf.loc]
            ## segments(x1, y1, f.range[2], delta.range[1], col = "red", lty = 2)
            segments(f.range[1], delta.range[2], f.range[2], delta.range[1], col = "red", lty = 2)
        }
    }
    points(cres$f[cres$icenter], cres$delta[cres$icenter], 
           col = mycols, pch = 19, cex = 1.1)
    if(centroids == "user")
        text(cres$f[cres$icenter], cres$delta[cres$icenter], labels = cres$icenter, pos = 1, cex = 0.6)
    ##--------------------
    ## Score plot
    ##--------------------
    if(plotType == "threeplots2"){
        plot(cres$hs, cres$scores, xlab = "h's", ylab = criterion, 
             main = "Selection with the best score")
        abline(v = cres$hs[cres$best.loc], col = "red", lty = 2)
    }
    
    ##--------------------
    ## Clustering results. First 2 principal comp if dim > 2.
    ##--------------------
    ndim <- ncol(dat)
    if(ndim > 2){
        pc2 <- prcomp(dat)$x[,c(1,2)]
        plot(pc2, col = mycols[cres$cluster],
             main = "First two principal components",
             xlab = "pc1", ylab = "pc2", pch = 1)
        points(pc2[cres$icenters,], col = mycols, pch = 4, cex = 3)
    }else{
        plot(dat, col = mycols[cres$cluster],
             main = "Original data",
             xlab = "x1", ylab = "x2", pch = 1)
        points(dat[cres$icenters,], col = mycols, pch = 4, cex = 3)
    }
    legend("topright",
           legend = c(paste0("cluster", 1:Final.nclust)),
           pch = 19,
           col = c(mycols[(1:Final.nclust - 1) %% length(mycols) + 1]),
           cex = 0.6
           )
    if(!is.null(file)) dev.off()
}
