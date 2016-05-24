`plotprintseq` <-
function(ngrid.r=4, ngrid.c=4,
           nspot.r=16, nspot.c=12, 
           gridorder=expand.grid(row=1:ngrid.c, col=1:ngrid.r),
           spotorder=list(x=nspot.r:1, y=nspot.c:1),
           rowmajor=FALSE, eps=1, delay1=100, delay2=2000){
    oldpar <- par(mar=par()$mar-c(0,0,2.5,0))
    on.exit(par(oldpar))
    plotpoints <- function(i, j, delay1=5000, delay2=10000){
      points(i+xy$x, j+xy$y, pch=15,
             cex=0.85, col="cyan")
      x <- 0
      for(k in 1:delay2)x <- x+1
      points(i+xy$x, j+xy$y, pch=15,
             cex=0.5, col="grey60")
      x <- 0
      for(k in 1:delay1)x <- x+1
    }

    xy <- gridorder-1
    names(xy) <- c("x","y")
    xy$x <- xy$x*(nspot.c+eps)
    xy$y <- xy$y*(nspot.r+eps)
    plot(c(1, ngrid.c*(nspot.c+eps)),
         c(1, ngrid.r*(nspot.r+eps)),
         type="n",xlab="",ylab="",axes=FALSE)
    mtext(side=1, line=1, 
          paste("Grid layout:  #rows of Grids =", ngrid.r,
                "   #columns of Grids =", ngrid.c))
    mtext(side=1, line=2.5, 
          paste("In each grid:  #rows of Spots =", nspot.r,
                "  #columns of Spots =", nspot.c))
    if (rowmajor)
      for(j in spotorder$x) for(i in spotorder$y)
        plotpoints(i,j, delay1=delay1, delay2=delay2)
    else
      for(i in spotorder$y) for(j in spotorder$x) 
        plotpoints(i,j, delay1=delay1, delay2=delay2)
  }

