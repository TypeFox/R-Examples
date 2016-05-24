
# functions for heatmap and color bars, thanks to Martin Sill

heat <- function(x, cols=maPalette(low="blue", mid="lightgrey", high="red", k=50), ylab, cex=1.5){
  # prevent outlier values to dominate the color scale -> restrict "maximal colors" to 2.5% and 97.5% quantiles
  qu <- quantile(x, c(0.025, 0.975), na.rm=TRUE)
  breaks <- seq(qu[1], qu[2], length.out=51)
  x[x < qu[1]] <- qu[1]
  x[x > qu[2]] <- qu[2]
  
  image(t(x)[1:ncol(x), nrow(x):1, drop=FALSE], col=cols, breaks=breaks, xlim=0.5 + c(0, ncol(x)),
        ylim = 0.5 + c(0, nrow(x)), x=c(1:ncol(x)), y=c(1:nrow(x)), axes=F, ylab="", xlab="", useRaster=TRUE)
  if(!missing(ylab))
    text(x=par("usr")[2], y=nrow(x):1, labels=ylab, xpd=T, cex=cex, pos=4)
}

fac2col <- function(x, cols=c("darkorange","darkred","thistle","cornflowerblue","olivedrab")){
  x <- as.factor(x)
  out <- rep("gray",length(x))
  lev <- levels(x)
  for(i in 1:length(lev)) 
    out[x==lev[i]] <- cols[i]
  out[is.na(x)] <- "white"
  names(out) <-  names(x)
  return(out)  
}

addfac <- function(x, ylab, cex=1.5){
  plot(y=c(0,ncol(x)), x=c(0,nrow(x)), type="n", ann=FALSE, xaxt="n", yaxt="n", xaxs="i", yaxs="i")
  for(i in ncol(x):1){
    ytop <- (ncol(x):1)[i]
    rect(xleft=seq(1, nrow(x))-1, xright=seq(1,nrow(x)), ytop=ytop, ybottom=ytop-1, border=NA, col=x[,i])
  }
  if(!missing(ylab))
    text(x=par("usr")[2], y=(ncol(x):1)-0.5, labels=ylab, xpd=T, cex=cex, pos=4)
}



## function to create a color legend matrix
legendmat <- function(data, Names, col.cont=maPalette(low="blue", mid="lightgrey", high="red", k=50),
                      lab.cex=1, col.ord=list(low="lightgreen", high="darkgreen"),
                      col.cat=c("darkorange","darkred","thistle","cornflowerblue","olivedrab","darkgrey","purple4","indianred","yellow2","darkseagreen4")){
  
  dc <- sapply(data, data.class)
  n.cont <- sum(dc == "numeric")
  n.ord <- sum(dc == "ordered")
  n.cat <- sum(dc == "factor")
  
  if(missing(Names))
    Names <- names(data)
  
  y <- max(n.ord, n.cat)
  K <- sapply(data[,dc == "factor"], function(x) length(levels(x)))
  mK <- 0
  if(n.cat != 0) 
    mK <- max(K)
  plot(0, xlim=c(-.5, 8 + mK), ylim=c(-.2, y), type="n", axes=FALSE, xlab="", ylab="")
  
  # legend for quantitative variables - just show 'min' and 'max' for any of them
  if(n.cont != 0){
    rect(.7, y-.7, 1.3, y+.2, col=col.cont[1])
    rect(1.7, y-.7, 2.3, y+.2, col=col.cont[length(col.cont)])
    text(.5, y-1, "quantitative", font=2, pos=2, xpd=TRUE, cex=lab.cex)
    text(c(1, 2), rep(y-1, 2), c("min","max"), cex=lab.cex)
  }
  
  # ordinal variables - also only show minimal and maximal categories (because color scheme is always developed within the low' and 'high' colors)
  if(n.ord != 0){
    rect(4.7, y-.7, 5.3, y+.2, col=col.ord$low)
    rect(5.7, y-.7, 6.3, y+.2, col=col.ord$high)
    
    N.cat <- names(data)[dc == "ordered"]    
    N.cat2 <- Names[dc == "ordered"]    
    for(i in 1:n.ord){
      L <- levels(data[,N.cat[i]])
      text(4.7, y-i, N.cat2[i], font=2, pos=2, cex=lab.cex)
      text(c(5,6), rep(y-i, 2), c(L[1],L[length(L)]), cex=lab.cex)
    }
  }
  
  # categorical variables
  if(n.cat != 0){
    N.cat <- names(data)[dc == "factor"]    
    N.cat2 <- Names[dc == "factor"]    
    for(i in 1:n.cat){
      rect(seq(8.7, length.out=mK), y-.7, seq(9.3, length.out=mK), y+.2, col=col.cat[1:mK])
      L <- levels(data[,N.cat[i]])
      l <- length(L)
      text(8.7, y-i, N.cat2[i], font=2, pos=2, cex=lab.cex)
      text(seq(9, length.out=l), rep(y-i, 2), L, cex=lab.cex)      
    }
  }
}
