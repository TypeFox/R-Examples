plot.emma <- function(x,n=50,fn,C=10,...)
{

  x1 <- seq(min(x$xspace[,1]), max(x$xspace[,1]), length = n)
  x2 <- seq(min(x$xspace[,2]), max(x$xspace[,2]), length = n)

  xx <- expand.grid(x1,x2)
  z <- matrix(fn(xx),nrow=length(x1))

  nrz <- nrow(z)
  ncz <- ncol(z)
  jet.colors <- colorRampPalette( c("red", "yellow") ) 
  nbcol <- 100
  color <- jet.colors(nbcol)
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  facetcol <- cut(zfacet, nbcol)

  for(i in 1:(C-1)){

    par(mfrow=c(1,2))

    res <- persp(x1, x2, z, col = color[facetcol], theta = 0, phi = 30, 
  	expand = 1, xlab = "x1", ylab = "x2", zlab = "f(x1, x2)", 
      ticktype = "detailed", main="Optimization function")
    if(i==1){ 
      ind1 <- 1
      ind2 <- x$nd
      points(trans3d(x$xpop[ind1:ind2,1], x$xpop[ind1:ind2,2], as.vector(x$ypop[ind1:ind2,1]), pmat = res), col = 1, bg="green", pch =21,cex=1)
      points(trans3d(x$xspace[x$Gb.arch[i],1], x$xspace[x$Gb.arch[i],2], as.vector(x$yspace[x$Gb.arch[i],1]), pmat = res), col = 1, bg="red", pch =21,cex=1)
    }

    res2 <- persp(x1, x2, z, col = color[facetcol], theta = 0, phi = 90, 
      expand = 1, xlab = "x1", ylab = "x2", zlab = "f(x1, x2)", 
      ticktype = "detailed", main=paste("Time instant t = ",i))
    if(i==1){
      ind1 <- 1
      ind2 <- x$nd
      points(trans3d(x$xpop[1:x$nd,1], x$xpop[1:x$nd,2], as.vector(x$ypop[1:x$nd,1]), pmat = res2), col = 1, bg="green", pch =21,cex=1)
      points(trans3d(x$xspace[x$Gb.arch[i],1], x$xspace[x$Gb.arch[i],2], as.vector(x$yspace[x$Gb.arch[i],1]), pmat = res2), col = 1, bg="red", pch = 21,cex=1)
    }
    if(i>1){ 
      ind1 <- ind2+1
      ind2 <- ind2+x$na
      points(trans3d(x$xpop[ind1:ind2,1], x$xpop[ind1:ind2,2], as.vector(x$ypop[ind1:ind2,1]), pmat = res2), col = 1, bg="green", pch =21,cex=1)
      points(trans3d(x$xspace[x$Gb.arch[i],1], x$xspace[x$Gb.arch[i],2], as.vector(x$yspace[x$Gb.arch[i],1]), pmat = res2), col = 1, bg="red", pch =21,cex=1)
    }
  }
}