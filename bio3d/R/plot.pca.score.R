"plot.pca.score" <-
function(x, inds=NULL, col=rainbow(nrow(x)), lab="", ... ) {

  # Produces a z-score plot for PC1 vs PC2,
  # PC3 vs PC2 and PC1 vs PC3 if given a
  # matrix "z"that contains the column wise
  # scores obtained from PCA 'pca.xyz'

  if(is.list(x))  x=x$z # output from pca.xyz()
  if(is.null(inds)) inds <- 1:nrow(x)
  
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  
  par(mfrow=c(2,2));par(pty="s")
  limits<-max(abs(c(range(x[,1]),range(x[,2]),range(x[,3]))))
  print(paste("axes limits: ",round(limits,2),sep=""))

  text.offset<-limits/19

  plot(x[,1], x[,2],  # pc1 vs pc2
       xlim=c(-limits,limits),
       ylim=c(-limits,limits),
       xlab="PC1", ylab="PC2",col=col, ...)  
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  text(x[inds,1]+text.offset,
       x[inds,2]+text.offset,
       labels = lab[inds])

  plot(x[,3], x[,2],  # pc3 vs pc3
       xlim=c(-limits,limits),
       ylim=c(-limits,limits),
       xlab="PC3", ylab="PC2", col=col, ...)
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  text(x[inds,3]+text.offset,
       x[inds,2]+text.offset,
       labels = lab[inds])

  plot(x[,1], x[,3],  # pc1 vs pc3
       xlim=c(-limits,limits),
       ylim=c(-limits,limits),
       xlab="PC1", ylab="PC3",col=col, ...)
  abline(h=0,col="gray",lty=2);abline(v=0,col="gray",lty=2)
  text(x[inds,1]+text.offset,
       x[inds,3]+text.offset,
       labels = lab[inds])

}
