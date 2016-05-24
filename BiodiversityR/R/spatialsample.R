`spatialsample` <-
function(x,method="random", n=5, xwidth=0.5, ywidth=0.5, xleft=0, ylower=0, xdist=0, ydist=0, plotit=T, plothull=F){
#    if (!require(splancs)) {stop("Requires package splancs")}
    xpos <- x[,1]
    ypos <- x[,2]
    minx <- min(xpos)
    maxx <- max(xpos)
    miny <- min(ypos)
    maxy <- max(ypos)
    xwidth <- xwidth/2
    ywidth <- ywidth/2
    if (method=="random") {    
        result <- array(dim=c(n,2))
        for (i in 1:n) {
            result[i,1] <- minx-1
            result[i,2] <- miny-1
            while((splancs::inout(cbind(result[i,1]-xwidth, result[i,2]-ywidth), x, bound=T)==F) || (splancs::inout(cbind(result[i,1]-xwidth, result[i,2]+ywidth), x, bound=T)==F) ||
                    (splancs::inout(cbind(result[i,1]+xwidth, result[i,2]-ywidth), x, bound=T)==F) || (splancs::inout(cbind(result[i,1]+xwidth, result[i,2]+ywidth), x, bound=T)==F)) {
                result[i,1] <- minx + (maxx-minx)*runif(1)
                result[i,2] <- miny + (maxy-miny)*runif(1)
            }
        }
    }
    if (method=="grid"  || method=="random grid" ) {
        if (xdist==0) {xdist <- (maxx-minx)/n}
        if (ydist==0) {ydist <- (maxy-miny)/n}
        if (xleft < minx) {xleft <- minx + xdist*runif(1)}
        if (ylower < miny) {ylower <- miny + ydist*runif(1)}
        a <- round((maxx-minx)/xdist)
        b <- round((maxy-miny)/ydist)
        result <- array(dim=c(a*b,2))
        for (i in 1:a) {
            for (j in 1:b) {
                result[((i-1)*b+j),1] <- xleft + (i-1)*xdist
                result[((i-1)*b+j),2] <- ylower + (j-1)*ydist
            }
        }
        i <- 1
        while (i <= nrow(result)) {
            if (splancs::inout(cbind(result[i,1]-xwidth, result[i,2]-ywidth), x, bound=T)==F) {
                result <- result[-i,]
            }else{
                i <- i+1
            }
        }
        i <- 1
        while (i <= nrow(result)) {
            if (splancs::inout(cbind(result[i,1]-xwidth, result[i,2]+ywidth), x, bound=T)==F) {
                result <- result[-i,]
            }else{
                i <- i+1
            }
        }
        i <- 1
        while (i <= nrow(result)) {
            if (splancs::inout(cbind(result[i,1]+xwidth, result[i,2]-ywidth), x, bound=T)==F) {
                result <- result[-i,]
            }else{
                i <- i+1
            }
        }
        i <- 1
        while (i <= nrow(result)) {
            if (splancs::inout(cbind(result[i,1]+xwidth, result[i,2]+ywidth), x, bound=T)==F) {
                result <- result[-i,]
            }else{
                i <- i+1
            }
        }
        if (n < nrow(result) && method=="random grid") {result <- result[(sample(nrow(result), n)),]}
    }
    if (plotit==T) {  
        graphics::rect(result[,1]-xwidth, result[,2]-ywidth, result[,1]+xwidth, result[,2]+ywidth)
        if (plothull==T) {
            points2 <- grDevices::chull(result)
            points3 <- c(points2, points2[1])
            graphics::lines(result[points3,])
        }
    }
    return(result)
}

