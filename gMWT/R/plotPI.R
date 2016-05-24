# Version: 30-11-2012, Daniel Fischer

plotPI <- function(X,g,type="pair",goi=NULL,mc=1,alg="Cnaive",col="black",highlight=NULL,hlCol="red",pch=20,zoom=FALSE,order=NULL,...){

  ifelse(is.null(order), chOrder <- TRUE, chOrder <- order)
  probs <- estPI(X=X,g=g,type=type,goi=goi,mc=mc,order=chOrder,alg=alg)

  plot(probs,col=col,highlight=highlight,hlCol=hlCol,pch=pch,zoom=zoom,...)
}