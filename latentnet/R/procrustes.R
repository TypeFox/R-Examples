.procrustes.Z.mean.C<-function(sample,Z.ref,center=FALSE,verbose=0){
  n<-dim(Z.ref)[1]
  G<-dim(sample[["Z.mean"]])[2]
  d<-dim(Z.ref)[2]
  S<-dim(sample[["Z"]])[1]
  ## Center Z.ref.
  Z.ref<-scale(Z.ref,scale=FALSE)

  Cret<-.C("procr_transform_wrapper",
           S=as.integer(S),
           n=as.integer(n),
           d=as.integer(d),
           G=as.integer(NVL(G,0)),
           Z.ref=as.double(Z.ref),
           Z=as.double(sample[["Z"]]),
           Z.mean=as.double(sample[["Z.mean"]]),
           verbose=as.integer(verbose),
           
           PACKAGE="latentnet")
  sample[["Z"]]<-if(d>0)array(Cret[["Z"]],dim=c(S,n,d))
  sample[["Z.mean"]]<-if(!is.null(G))array(Cret[["Z.mean"]],dim=c(S,G,d))
  
  sample
}

.procr <- function(x, ...) UseMethod(".procr")
.procr.matrix <- function(x, ref, ..., scale=FALSE, reflect=TRUE){
  ref <- sweep(ref, 2, colMeans(ref), "-")
  x <- sweep(x, 2, colMeans(x), "-")

  M <- crossprod(x, ref)
  M.svd <- svd(M)
  R <- (if(reflect) M.svd$u%*%t(M.svd$v) else M.svd$u%*%diag(c(det(M.svd$u%*%t(M.svd$v)),rep(1,ncol(ref)-1)),nrow=ncol(ref))%*%t(M.svd$v)) * if(scale) sqrt(sum(ref^2)/sum(x^2)) else 1
  R
}

.procr.ergmm.model <- function(x, A, ref, ...){
  .procr(A, ref, scale="scaling" %in% latent.effect.invariances[[x[["latentID"]]]],reflect="reflection" %in% latent.effect.invariances[[x[["latentID"]]]])
}
