plot.dmap <- function(x,...){

  if(x$neigen>=3){
    scatterplot3d(x$X[,1],x$X[,2],x$X[,3],pch=20,cex.symbols=.6,xlab="1st Diffusion Coord.",ylab="2nd Diffusion Coord.",zlab="3rd Diffusion Coord.",main="3-Dimensional Diffusion Map")
  }else if(x$neigen==2){
    plot(x$X[,1],x$X[,2],pch=20,cex=.6,xlab="1st Diffusion Coord.",ylab="2nd Diffusion Coord.",main="2-Dimensional Diffusion Map")
  }else {
    stripchart(x$X,method='jitter',pch=20,cex=.6,xlab="1st Diffusion Coord.",main="1-Dimensional Diffusion Map")
  }
  
}
