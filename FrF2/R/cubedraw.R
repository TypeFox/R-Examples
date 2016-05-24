cubedraw <- function(min = 0, max = 1, cex = 4, cex.lab=2, cex.axlab=2, center=FALSE,star=FALSE,a=1,pch=20,
     pchextra=pch, main=NULL, sub=NULL, ticks=FALSE, xlab="A", ylab="B", zlab="C",
     xends=c(-1,1),yends=c(-1,1),zends=c(-1,1), 
     x.ticklabs=c(xends[1],NA,xends[2]), y.ticklabs=c(yends[1],NA,yends[2]),
     z.ticklabs=c(zends[1],NA,zends[2]), zadj=c(0.8,0.9), nint=c(2,2,2),...)
 {
 ## Purpose:
 ## Draw cube with corners and potentially center and or stars
 ## flexible regarding plot symbols
 ## and colors
    oldpar <- par()
    par(mar=c(5, 3, 4, 3) + 0.1, xpd=TRUE)
    cub <- expand.grid(A=c(min,max),B= c(min,max),C= c(min,max))
    res3d<-myscatterplot3d(cub,cex.symbols=cex,pch=pch,grid=FALSE,box=FALSE,
         main=main, sub=sub, tick.marks=ticks, xlab=xlab,ylab=ylab,zlab=zlab,
         cex.lab=cex.lab,cex.axis=cex.axlab, xlim=c(-1,1), ylim=c(-1,1), zlim=c(-1,1),
         x.ticklabs=x.ticklabs, y.ticklabs=y.ticklabs,
         z.ticklabs=z.ticklabs, zadj=zadj, lab=nint,...)
    mid <- (max+min)/2
    ## hidden corner + lines
    res3d$points3d(cub[c(1,3,4,3,7), ], cex = cex, type = 'b', lty = 1,lwd=3,col=grey(0.4),pch=pch)
    ## visible corners + lines:
    res3d$points3d(cub[c(2,1,5,6,2,4,8,6,5,7,8) ,], cex = cex, type = 'b', lty = 1,lwd=3,pch=pch)
    res3d$points3d(cub,cex=cex,pch=pch,lwd=3) ## points in top layer
    if(center){
        xnullvis<-rbind(c(mid,min,min),c(mid,min,max),c(mid,max,max))
        xnullinvis <- rbind(c(mid,max,max),c(mid,max,min),c(mid,min,min))
        ynullvis<-rbind(c(max,mid,min),c(max,mid,max),c(min,mid,max))
        ynullinvis <- rbind(c(min,mid,max),c(min,mid,min),c(max,mid,min))
        znullvis<-rbind(c(min,min,mid),c(max,min,mid),c(max,max,mid))
        znullinvis <- rbind(c(max,max,mid),c(min,max,mid),c(min,min,mid))
        res3d$points3d(rbind(c(mid,mid,mid)), cex = cex, pch=pchextra, type = 'p')
        res3d$points3d(xnullvis, cex = cex, type = 'l', lty = max,lwd=max)
        res3d$points3d(ynullvis, cex = cex, type = 'l', lty = max,lwd=max)
        res3d$points3d(znullvis, cex = cex, type = 'l', lty = max,lwd=max)
        res3d$points3d(xnullinvis, cex = cex, type = 'l', lty = 3,lwd=max)
        res3d$points3d(ynullinvis, cex = cex, type = 'l', lty = 3,lwd=max)
        res3d$points3d(znullinvis, cex = cex, type = 'l', lty = 3,lwd=max)
        xc<-rbind(c(min,mid,mid),c(mid,mid,mid),c(max,mid,mid))
        yc<-rbind(c(mid,min,mid),c(mid,mid,mid),c(mid,max,mid))
        zc<-rbind(c(mid,mid,min),c(mid,mid,mid),c(mid,mid,max))
        res3d$points3d(xc, cex = cex, type = 'l', lty = 3,lwd=max)
        res3d$points3d(yc, cex = cex, type = 'l', lty = 3,lwd=max)
        res3d$points3d(zc, cex = cex, type = 'l', lty = 3,lwd=max)
    }
    if (star)
    {
        xc<-rbind(c(-a,mid,mid),c(mid,mid,mid),c(a,mid,mid))
        yc<-rbind(c(mid,-a,mid),c(mid,mid,mid),c(mid,a,mid))
        zc<-rbind(c(mid,mid,-a),c(mid,mid,mid),c(mid,mid,a))
        res3d$points3d(xc, cex = cex, pch=pchextra, type = 'b', lty = 3,lwd=3)
        res3d$points3d(yc, cex = cex, pch=pchextra, type = 'b', lty = 3,lwd=3)
        res3d$points3d(zc, cex = cex, pch=pchextra, type = 'b', lty = 3,lwd=3)
    }
    par <- oldpar
    list(res3d=res3d, cub=cub)
 }