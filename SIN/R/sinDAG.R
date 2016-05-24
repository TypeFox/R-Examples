## order has to be a list, e.g. as.list(1:p)
sinDAG <- function(order, S, n, holm=TRUE){
  order <- unlist(order)
  old.order <- dimnames(S)[[1]]
  S <- S[order,order]
  p <- dim(S)[1]
  pvals <- diag(rep(1,p))
  labels <- dimnames(S)[[1]]
  dimnames(pvals) <- list(labels, labels)
  for(i in 2:p){
    corr.part <- -cov2cor(solve(S[1:i,1:i]))
    pvals[i,1:(i-1)] <- pvals[1:(i-1),i] <-
      simpvalueVec(c(corr.part[i,1:(i-1)]),n-i-1,p)
  }
  for(i in 1:p) pvals[i,i] <- NA
  if(holm==TRUE) pvals <- holm(pvals)
  return(zapsmall(pvals[old.order,old.order]))
}
plotDAGpvalues <- function(pvals, legend=TRUE, legendpos=NULL){
  createDAGlabels <- function(pvals){
    p <- dim(pvals)[1]
    labels <- c()
    for(i in 2:p){
      for(j in 1:(i-1)){
        labels <- c(labels, paste(j,i, sep="->"))
      }
    }
    return(labels)
  }
  par(mar=c(6,5,2,2)+0.1)
  DAGlab <- createDAGlabels(pvals)
  DAGpvals <- pvals[upper.tri(pvals)] 
  temp <- length(DAGlab)
  plot(as.factor(1:temp), DAGpvals, type="n",
       ylab="P-value", xlab="", axes=FALSE, ylim=c(0,1), cex.lab=1.2, las=2)
  title(xlab = "Edge", line = 4, cex.lab=1.2)
  axis(1, at=1:temp, labels=DAGlab[1:temp], las=2)
  axis(2, at=seq(0,1,by=0.1), las=1)
  temp2 <- sapply(seq(0,1,by=0.1), abline, 0, lty="dotted", col="grey")
  plot(as.factor(1:temp), DAGpvals, add=TRUE, axes=FALSE)
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
