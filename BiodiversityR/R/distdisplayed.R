`distdisplayed` <-
function(x,ordiplot,distx="bray",plotit=T,addit=F,method="spearman",permutations=100,abline=F,gam=T,...) {
    if (gam==T) {
#        if (!require(mgcv)) {stop("Requires package mgcv")}
    }
    if (inherits(x, "dist")) {
        dist1 <- x
        xlab <- attr(x,"method")
        if (is.null(xlab)) {
            xlab <- "distance matrix 1"
        }else{
            xlab <- paste(xlab, "distance")
        }
    }else{
        dist1 <- vegdist(x, method = distx)
        xlab <- "original distance"
    }
    if (inherits(ordiplot, "dist")) {
        dist2 <- ordiplot
        ylab <- attr(dist2,"method")
        if (is.null(ylab)) {
            ylab <- "distance matrix 2"
        }else{
            ylab <- paste(ylab, "distance")
        }
    }else{
        ordiscores <- scores(ordiplot,display="sites")
        dist2 <- vegdist(ordiscores,method="euclidean")
        ylab <- "distance in ordination plot"
    }
    if (plotit==T) {
        if (addit==F) {
            graphics::plot(dist1, dist2, xlab=xlab, ylab=ylab)
        }
        if (abline==T) {abline(0,1)}
        if (gam==T){
            data <- data.frame(cbind(dist1,dist2))
            seq <- order(data[,1])
            sorted <- data
            sorted[1:nrow(data),] <- data[seq,]
            gamresult <- gam(dist2~s(dist1),data=sorted)
            newdata <- data.frame(seq(min(sorted[,1]), max(sorted[,1]), length = 1000))
            colnames(newdata) <- "dist1"
            gamresult2 <- predict(gamresult,newdata) 
            graphics::points(newdata$dist1,gamresult2,type="l",lwd=2,col="red")
        }
    }
    result2 <- mantel(dist1,dist2,method=method,permutations=permutations,...)
    if (gam==T) {
        return(list(gamanalysis=summary(gamresult),mantelanalysis=result2))
    }else{
        return(result2)
    }
}

