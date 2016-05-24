"plot.voronoi.polygons" <- function(x,which,color=TRUE,...){
  lx <- length(x)
  if(missing(which))
    which <- 1:lx
  lmax <- function(x)
    apply(x,2,max)
  lmin <- function(x)
    apply(x,2,min)
  lmean <- function(x)
    apply(x,2,mean)
  xy.max <- apply(sapply(x,lmax),1,max)
  xy.min <- apply(sapply(x,lmin),1,min)
  xy.mean <- sapply(x,lmean)
  plot(x[[which[1]]],type="n",xlim=c(xy.min["x"],xy.max["x"]),
       ylim=c(xy.min["y"],xy.max["y"]),...)
  colors <- heat.colors(lx)
  for(i in which){
    polygon(x[[i]],col=colors[i])
    text(xy.mean[,i]["x"],xy.mean[,i]["y"],i)
  }
  title(paste("plot of",deparse(substitute(x))))
}
