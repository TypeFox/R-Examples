plotjointICC.dRm <- function(
  object,
  item.subset = "all",
  legend = TRUE,
  xlim = c(-4, 4),
  ylim = c(0, 1),
  xlab = "Latent Dimension",
  ylab = "Probability to Solve",
  lty = 1,
  legpos="topleft",
  main = "ICC plot",
  col = NULL,
  ...
){
# produces one common ICC plot for Rasch models only
# object of class "dRm"
# item.subset...specify items that have to be plotted; if NA, all items are used
# legend...if legend should be plotted

  theta <- seq(xlim[1L], xlim[2L], length.out = 201L)

  if(any(item.subset=="all")){
    it.legend <- 1:dim(object$X)[2]
  } else {
    if(is.character(item.subset)){
      it.legend <- item.subset
      betatemp  <- t(as.matrix(object$betapar))
      colnames(betatemp) <- colnames(object$X)
      object$betapar <- betatemp[,item.subset]
    } else {
      it.legend <- colnames(object$X)[item.subset]
      object$betapar <- object$betapar[item.subset]
    }
    object$X <- object$X[,item.subset]                            #pick out items defined in itemvec
  }

  th.ord <- order(theta)

  p.list <- plist.internal(object, theta)
  p.list <- lapply(p.list, function(x){ x[,-1L] })               #Delete 0-probabilites
  p.mat  <- matrix(unlist(p.list), ncol = length(p.list))
  text.ylab <- p.mat[(1:length(theta))[theta==median(theta)],]

  if(is.null(main)) main=""
  if(is.null(col)) col=1:(dim(p.mat)[2])
  #pmICCs<-cbind(sort(theta),p.mat[th.ord,])
  matplot(sort(theta),p.mat[th.ord,],type="l",lty=lty,col=col,
          main=main,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
  if(length(object$betapar)>20) old_par <- par(cex=0.7) else old_par <- par(cex=1)
  on.exit(par(old_par))
  
  if(is.character(legpos)){
    if(!legend){
      sq <- seq(0.65,0.35,length.out=length(object$betapar))
      x  <- qlogis(sq,sort(-object$betapar))
      text(x=x,y=sq,labels=it.legend[order(-object$betapar)],col=col[order(-object$betapar)],...)
    } else {
      legend(legpos,legend=paste("Item",it.legend[order(-object$betapar)]),lty=lty, col=col[order(-object$betapar)],...)
    }
  }
}
