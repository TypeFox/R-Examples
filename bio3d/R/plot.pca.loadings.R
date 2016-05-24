"plot.pca.loadings" <-
function(x, resnums= seq(1,(length(x[,1])/3), 25), ... ) {

  # Plot residue loadings along PC1 to PC3 if given an xyz
  # C-alpha matrix of "loadings" (e.g. as returned from
  # 'pca.xyz' such a 'pca.trj$loadings')
  # For more info see 'pca.res.loadings'
  #
  # To Do: add gap.cols options

  if(is.list(x))  x=x$U
  
  pos <- resnums*3
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  par(mfrow=c(3,1), mar=c(4,4,2,2))

  plot(abs(x[,1]),main="",type="h",
       axes=FALSE, xlab="Index Number", ylab="PC1")
  axis(1, at=pos, labels=resnums)
  axis(2)
  box()

  plot(abs(x[,2]),main="",type="h",
       axes=FALSE,xlab="Index Number",ylab="PC2")
  axis(1, at=pos,labels=(pos)/3)
  axis(2)
  box()

  plot(abs(x[,3]),main="",type="h",
       axes=FALSE,xlab="Index Number",ylab="PC3")
  axis(1, at=pos,labels=(pos)/3)
  axis(2)
  box()
}
