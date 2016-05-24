slidevector <-
function(data,dim=2,itmax=125,eps=1e-12){
  if (sum(data<0) > 0) stop("data for the slide-vector model should be positive")
#check for number of rows equals number of columns
  if (nrow(data)!=ncol(data)) stop("the same number of rows and columns are expected")

    hjw <- .jnew("asymmetry/slidevector")     # create instance of slidevector class
  dis <- .jcall(hjw, "[[D","slidevectormodel",.jarray(data,dispatch=TRUE),as.integer(dim),as.integer(itmax),as.double(eps))
  ndim <- .jcall(hjw,"I","getDimension")
  nobs <- .jcall(hjw,"I","getNobs")
  niter <- .jcall(hjw,"I","getNiter")
  stress <- .jcall(hjw,"D","getStress")
  X=t(sapply(dis,.jevalArray))
  mat<-X[1:nobs,]
  Z<-X[nobs+1,]
  if(!is.null(rownames(data)))
    rownames(mat)<-rownames(data)
  colnames(mat) <- paste("D",1:(dim(mat)[2]),sep="")

  result<-list(ndim=ndim,stress=stress,confi=mat,slvec=Z,niter=niter,nobs=nobs)
  class(result)<-"slidevector"
  return(result)
}
