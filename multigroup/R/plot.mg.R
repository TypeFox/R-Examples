#' @title Plots for multigroup objects 
#' @description plots of variables (loadings) and individuals (scores) if TRUE
#' @param x results of multigroup method in the package
#' @param axes by default the first two components
#' @param cex character expansion for text by default .85
#' @param font.lab type of font by default 3
#' @param \dots Further arguments are ignored
#' @return loadings and scores plots
#' @S3method plot mg
#' @method plot mg
plot.mg <- function(x, axes=c(1,2), cex=NULL, font.lab= NULL, ...){
  #=========================================================================
  #                            Preparing inputs
  #=========================================================================
  Group=x$Group
  xax=axes[1]
  yax=axes[2]
  lab.x <- paste("Dim ", axes[1],sep = "")
  lab.y <- paste("Dim ", axes[2],sep = "")
  if(is.null(cex)) {cex=.85}
  if(is.null(font.lab)) {font.lab=3}
   cooll <- function(n, alfa=1) {
    hues = seq(15,375,length=n+1)
    hcl(h=hues, l=65, c=100, alpha=alfa)[1:n]
  }
  #=========================================================================
  #                         Common loadings plot
  #=========================================================================
  AA=x$loadings.common
  PQ=nrow(AA)
  w1=AA[,xax]
  w1=matrix(w1, ncol=1)
  w2=AA[,yax]
  w2=matrix(w2, ncol=1)
  vv=c(rep(0,PQ))
  uu=c(rep(0,PQ))
  minlimx   <- min(c(w1,w2))
  maxlimx   <- max(c(w1,w2))
  lim =c(minlimx,maxlimx)
  plot(w1,w2, type="n",ylim=lim ,xlim=lim ,xlab =lab.x, ylab=lab.y,main="Loadings plot",asp= 1)
  abline(h = 0, v = 0, col= "gray60")
  www=cbind(w1,w2)
  text(www, labels=rownames(AA), cex=cex, font.lab= font.lab) 
  #=========================================================================
  #                            score plot
  #=========================================================================
  dev.new()
  TT=x$Con.Data %*% x$loadings.common
  TT=TT[,c(1,2)]
  cp1 <- round(var(TT[,xax])/sum(diag(var(TT))), digits = 2)
  cp2 <- round(var(TT[,yax])/sum(diag(var(TT))), digits = 2)
  lab.x <- paste("Dim ", xax, sep = "")
  lab.y <- paste("Dim ", yax, sep = "")
  labs_col = cooll(nlevels(Group), alfa=1)
  rep.labs_col=rep(labs_col,as.vector(table(Group)))
  plot(TT[,xax],TT[,yax], xlab = lab.x, ylab = lab.y,type="n",main="Individual plot") 
  abline(h = 0, v = 0, col= "gray60")
  text(TT,labels=rownames(TT), cex=cex, font.lab=font.lab, col=rep.labs_col)
}
