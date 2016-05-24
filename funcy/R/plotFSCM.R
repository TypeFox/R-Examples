#
#  Copyright (C) 2011-2015 Christina Yassouridis
#  
                                        #

plotOverview <- function(object, showLegend){
    
    temp <- object@trends$temp
    k <- max(object@cluster)
    time <- object@time
    
    ctrs1 <- temp[,-1, drop=FALSE]
    ctrs2 <- object@centers
    
    dist <- matrix(NA,k,k)
    for(i in 1:k)
        for(j in 1:k)
            dist[i,j] <- sum((ctrs1[,i]-ctrs2[,j])^2)
    labels <- resort(dist)

    reset <- function() {
        par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
        plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
    }

    .pardefault <- par() 
    par(xpd=T, mfrow=c(2,2) , mar=c(3,1,3,10))
    
    plotClusterLocs(object, showLegend=FALSE)
    

    matplot(time, ctrs1, type='l', lty=1, col=labels[2,],  ylab="",ylim=c(min(temp[,1],ctrs1),max(temp[,1],ctrs1)))
    title("Temporal cluster trends",line=1)
    lines(time, temp[,1],type='l', lty=1, col=k+1, lwd=2)
    matplot(time, object@centers, type='l', lty=1, col=labels[1,])
    title("Temp + Overall", line=1)

    plotSpatials(object, showLegend=showLegend)

    if(showLegend){
        reset()
        legend("top", c("overall trend", paste("cluster",1:k)), col =
                   c(k+1,labels[1,]),pch="o", bty = "n")
    }
    
    
}

plotClusterLocs <- function(object, showLegend){
    k <- max(object@cluster)
    location <- object@location
    cl <- as.factor(object@cluster)
    levels(cl) <- col <- 1:k
    cl <- as.character(cl)
    plot(location[,1], location[,2], col=cl, type="p", pch=16, ylab="")
    title("Cluster location", line=1)
    if(showLegend)
        legend("topright", inset=c(-0.5,0), paste("cluster", 1:k), col =
              col, pch="o", bty = "n")
}


plotSpatials <- function(object, showLegend){
    location <- object@location
    spatial <- object@trends$spatial
    if(!is.null(dim(spatial)))
            spatial <- spatial[1,]

    N <- dim(location)[1]
      
    a <- (spatial+max(abs(spatial)))
    spatial2col <- as.numeric((a-min(a))/(max(a)-min(a)))
    rbPal <- colorRampPalette(c('yellow','blue'))
    Col <- rbPal(N)[as.numeric(cut(spatial2col,breaks = N))]
    legend_image <- as.raster(matrix(rbPal(20), ncol=1))

    ##layout(matrix(1:2,ncol=2), width = c(2,1),height = c(1,1))
    plot(location[,1],location[,2],col=Col,cex=1.5,pch=20, xlab="x", ylab="y")
    title(main="Spatial dependency", line=1)
    mtext(paste("Method", object@methodName), outer=TRUE)

    if(showLegend){
        x <- location[,1]
        y <- location[,2]
        relx <- ceiling(diff(range(location[,1]))/20)
        rely <- ceiling(diff(range(location[,2]))/5)
        rasterImage(legend_image, max(x)-relx, max(y)-rely,max(x),max(y))
        
        text(x=rep(max(x)-0.5*relx,5), y=seq(max(y)-rely,max(y),l=5), labels =
                 seq(1,0,l=5), cex=0.7, lwd=2)
    }
}


plotDeviations <- function(object, showLegend){
    location <- object@location
    N <- dim(location)[1]
    devScaled <- object@prms$sigma
    rbPal <- colorRampPalette(c('lightblue','darkblue'))
    Col <- rbPal(N)[as.numeric(cut(devScaled, breaks = N))]
    eqscplot(location[,1],location[,2],pch = 20, cex=1.5, col = Col,main="Deviations")
}


plotLoc <- function(object, showLegend){
    location <- object@location
    nc <- dim(location)[1]
    cex <- ifelse(nc>100, 100/nc, 0.7)
    par(mfrow=c(1,2))
    eqscplot(location[,1],location[,2], xlab="location 1",
             ylab="location 2", cex=cex)
    textxy(location[,1],location[,2],1:length(location[,1]), cex=max(cex,0.5))
    
    B=stationary.cov(location, location, Covariance="Matern")
    B <- as.dist(B)
    check <- cmdscale(1-B,k=2)
    sunflowerplot(check[,1],check[,2],
             main="Matern covariance distance plot",
             ylab="md-scaled variable 2", xlab="md-scaled variable 1",
                  cex=cex)
    textxy(check[,1], check[,2], 1:length(check[,1]), cex=max(cex,0.5))
    par(mfrow=c(1,1))

    ##Get the order of the closest points and verify with plot
    B <- as.matrix(B)
    a <- matrix(0,nc,nc)
    for(j in 1:nc)
        for(i in 1:nc) a[i,j]=which(B[,j]==sort(B[,j],
                                                decreasing=T)[i])[1]
    return(a)

}

