sinUG <- function(S,n, holm=TRUE){
  p <- dim(S)[1]
  pvals <- simpvalueMx(-cov2cor(solve(S)), n-p-1, p)
  dimnames(pvals) <- dimnames(S)
  if(holm==TRUE) pvals <- holm(pvals)
  return(zapsmall(pvals))
}
plotUGpvalues <- function(pvals, legend=TRUE, legendpos=NULL){
#   vecUGpvalues <- function(pvals){
#     p <- dim(pvals)[1]
#     pvec <- c()
#     for(i in 1:(p-1)){
#       pvec <- c(pvec, pvals[i,(i+1):p])
#     }
#     return(pvec)
#   }
  createUGlabels <- function(pvals){
    p <- dim(pvals)[1]
    labels <- c()
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        labels <- c(labels, paste(i,j, sep="-"))
      }
    }
    return(labels)
  }
  par(mar=c(6,5,2,6)+0.1)
  UGlab <- createUGlabels(pvals)
  UGpvals <- pvals[lower.tri(pvals)] #vecUGpvalues(pvals)
  temp <- length(UGlab)
  plot(as.factor(1:temp), UGpvals, type="n",
       ylab="P-value", xlab="", axes=FALSE, ylim=c(0,1), cex.lab=1.2, las=2)
  title(xlab = "Edge", line = 4, cex.lab=1.2)
  axis(1, at=1:temp, labels=UGlab[1:temp], las=2)
  axis(2, at=seq(0,1,by=0.1), las=1)
  temp2 <- sapply(seq(0,1,by=0.1), abline, 0, lty="dotted", col="grey")
  plot(as.factor(1:temp), UGpvals, add=TRUE, axes=FALSE)
  box()
  p <- dim(pvals)[1]
  plotlabels <- as.character(1:p)
  for(i in 1:p){
    plotlabels[i] <- paste(plotlabels[i],dimnames(pvals)[[1]][i], sep="  ")
  }
  if(legend==TRUE){
    if(is.null(legendpos)){
      legend(temp+0.6,1, x.intersp=-0.3, 
             plotlabels, bg="white", xjust=1, yjust=1)
    }
    else{
      x <- legendpos[1]
      y <- legendpos[2]
      legend(x,y, x.intersp=-0.3, 
             plotlabels, bg="white", xjust=1, yjust=1)
    }
  }
}
