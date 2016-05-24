plotObserved <- function(data, results, pch=NULL, col=NULL, xlab=NULL, ylab=NULL, main=NULL, SE=T){
  genome.perc <- read.csv(file = system.file("genome.scale.csv", package = "SAGA"))[,2:3]
  x <- cbind(data, genome.perc[data[,1],2])
  colnames(x)[4] <- "gen.perc"
  if(is.null(pch)) pch <- 16
  if(is.null(col)) col <- "black"
  if(is.null(xlab)) xlab <- "% P1 Genome"
  if(is.null(ylab)) ylab <- "Phenotype Measure"
  if(length(unique(x[,4])) != length(x[,4])){
    xvals <- jitter(x[,4])
  }else{
    xvals <- x[,4]
  }
  if(SE==T){
    yvals <- vector()
    for(i in 1:nrow(x)){
      high <- sum(x[i,2:3])
      if(i == 1) yvals[i] <- sum(x[i,2:3])
      low <- x[i,2] - x[i,3]
      if(high > max(yvals)) yvals <- c(yvals, high)
      if(low  < min(yvals)) yvals <- c(yvals, low)
    }
    high <- max(yvals)
    low <- min(yvals)
    plot(x=xvals, y=x[,2], ylab=ylab, xlab=xlab, xaxt="n", pch=pch, main=main, ylim=c(low,high))
    for(i in 1:nrow(x)){
      lines(x=rep(xvals[i], 2), y= c(sum(x[i,2:3]), x[i,2] - x[i,3]))
    }
  }else{
    plot(x=xvals, y=x[,2], ylab=ylab, xlab=xlab, xaxt="n", pch=pch, main=main)
  }
  axis(side=1,labels=c(0,50,100), at=c(0,50,100))
  abline(glm(x[,2]~x[,4], weights = data[, 2] ^ - 2), lty="dashed", col="blue")
}