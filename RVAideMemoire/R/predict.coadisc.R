# ade4 : as.dudi, dist.quant

predict.coadisc <-
function(object,newdata,dim=object$nf,method=c("mahalanobis", "euclidian"),...) {
  if (dim>object$nf) {stop("'dim' > number of dimensions in CDA")}
  method <- match.arg(method)
  df <- eval.parent(object$call$df)
  if (ncol(df)!=ncol(newdata)) {stop("incorrect dimensions of 'newdata'")}
  N <- sum(df)
  df <- df/N
  col.w <- colSums(df)
  fa <- as.matrix(object$fa)[,1:dim]
  coord <- matrix(0,ncol=dim,nrow=nrow(newdata))
  colnames(coord) <- paste0("DS",1:dim)
  newdata <- as.matrix(newdata)
  for (i in 1:nrow(newdata)) {
    new <- t(newdata[i,])/N
    row.w.new <- apply(new,1,sum)
    new <- new/row.w.new
    new <- sweep(new,2,col.w)
    X <- ade4::as.dudi(as.data.frame(new),1/col.w,row.w.new,scannf=FALSE,nf=dim, 
	call=match.call(),type="coarp",full=TRUE)
    tab <- as.matrix(X$tab)
    coord[i,] <- tab%*%fa
  }
  Y <- eval.parent(object$call$fac)
  li <- if (dim==1) {t(t(object$li[,1:dim]))} else {object$li[,1:dim]}
  centr <- apply(li,2,function(x) tapply(x,Y,mean))
  res <- character(nrow(coord))
  for (i in 1:nrow(coord)) {
    r <- rbind(coord[i,],centr)
    dis <- if (method=="mahalanobis") {
	as.matrix(ade4::dist.quant(r,method=3))[-1,]
    } else {
	as.matrix(dist(r,method="euclidian"))[-1,]
    }
    res[i] <- rownames(dis)[which.min(dis[,1])]
  }
  res <- factor(res)
  return(res)
}
