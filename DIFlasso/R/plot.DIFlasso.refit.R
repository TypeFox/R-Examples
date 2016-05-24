plot.DIFlasso.refit <-
function(x, decreasing = TRUE, ...){
  oldpar <- par(no.readonly=TRUE)
  plot.mat <- x$gamma
  dif.strength <- rowSums(plot.mat^2)
  if(decreasing){
    plot.mat <- plot.mat[order(dif.strength,decreasing = decreasing), ,drop = FALSE]
  }
  y.min <- min(plot.mat)
  y.max <- max(plot.mat)
  par(mar=oldpar$mar+c(0,0,0,2))
  plot(plot.mat[,1],ylim=c(y.min,y.max),type="b",xaxt="n",
       ylab="",xlab="covariates", main="Item-specific parameter estimates", ...)
  for(i in 1:ncol(plot.mat)){
    lines(plot.mat[,i],type="b")
  }
  abline(h=0,lty=2)
  axis(1,labels=dimnames(plot.mat)[[1]],at=1:nrow(plot.mat))
  axis(4,labels=dimnames(plot.mat)[[2]],at=plot.mat[nrow(plot.mat),],las=2)
  on.exit(par(oldpar))
}
