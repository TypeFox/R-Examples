impute=function(object,i,j,unscale=TRUE){
  ###obj is a svd object produced by our Softimpute algorithms
  v=as.matrix(object$v)
 vd=v*outer(rep(1,nrow(v)),object$d)
 out= suv(as.matrix(object$u),vd,i,j)
 if(unscale){
   biats=attributes(object)
   if(any(grep("biScale",names(biats)))){
     out=out*biats[["biScale:row"]]$scale[i]*biats[["biScale:column"]]$scale[j]
     out=out+biats[["biScale:row"]]$center[i]+biats[["biScale:column"]]$center[j]
   }
 }
 out
}
