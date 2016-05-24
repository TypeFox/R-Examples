#' @title loadings plot of X and Y
#' @description plots of variables (loadings)
#' @param X common loadingsassociated with X
#' @param Y common loadingsassociated with Y
#' @param INERTIE if there is information about inertia
#' @param axes a vector of two selected components 
#' @param cex character expansion for text by default .85
#' @param font.lab type of font by default 3
#' @return loadings plot
#' @export
#' 
#' @examples
#' data(oliveoil)
#' DataX = oliveoil[,2:6]
#' DataY = oliveoil[,7:12]
#' Group = as.factor(oliveoil[,1])
#' res.mgPLS = mgPLS (DataX, DataY, Group)
#' X=res.mgPLS$loadings.commo$X; Y=res.mgPLS$loadings.commo$Y
#' loadingsplotXY(X, Y, axes=c(1,2), INERTIE=res.mgPLS$noncumper.inertiglobal)
loadingsplotXY <- function(X, Y,  axes=c(1,2), INERTIE=NULL, cex=NULL, font.lab= NULL){
  
  
  if (max(axes)>ncol(X))
    stop("\nOops one of the axes value is greater than number of existing dimensions")
  
  if (max(axes)>ncol(Y))
    stop("\nOops one of the axes value is greater than number of existing dimensions")
  
  
  if(is.null(rownames(X))) {
    rownames(X) = paste('X', 1:nrow(X), sep='')
  }
  
  
  if(is.null(rownames(Y))) {
    rownames(Y) = paste('Y', 1:nrow(Y), sep='')
  }
  
  
  AA = rbind(X[,axes],Y[,axes])
  P  = nrow(X)
  Q  = nrow(Y)
  PQ = P+Q
  w1 = AA[,1]
  w1 = matrix(w1, ncol=1)
  w2 = AA[,2]
  w2 = matrix(w2, ncol=1)
  vv = c(rep(0,PQ))
  uu = c(rep(0,PQ))
  
  minlimx   <- min(c(w1,w2))
  maxlimx   <- max(c(w1,w2))
  lim =c(minlimx   ,maxlimx   )
  
  #-----------------------
  xax=axes[1]
  yax=axes[2]
  
  
  Dim11 = paste("Dim ", axes[1], sep = "")
  Dim21 = paste("Dim ", axes[2], sep = "")
  
  
  Dim12 = paste(Dim11, INERTIE[axes[1]], sep =" (")
  Dim22 = paste(Dim21, INERTIE[axes[2]], sep =" (")
  
  lab.x = paste(Dim12, "%)", sep="")
  lab.y = paste(Dim22, "%)", sep="")
  #------------------------
  
  plot(w1, w2, type="n", ylim=lim ,xlim=lim ,xlab = lab.x, ylab = lab.y)
  abline(h = 0, v = 0, , col= "gray60")
  arrows(vv, uu, w1, w2,lwd=2,length=.2, lty=c(rep(1,P)), col=rep(c(1,2),c(P,Q)))
  www = cbind(w1,w2)
  text(www,labels=rownames(AA), cex=1, font.lab= 3, pos=4, , col=rep(c(1,2),c(P,Q))) 
}
