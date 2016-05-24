nystrom <- function(dmap,Dnew,sigma=dmap$epsilon){
  
  Nnew = dim(Dnew)[1]
  Nold = dim(Dnew)[2]

  if(Nold != dim(dmap$X)[1]){
    print("dimensions don't match")
    break
  }

  Xnew = exp(-Dnew^2/(sigma))
  v = apply(Xnew, 1, sum)
  Xnew = Xnew/matrix(v,Nnew ,Nold)
  #nystrom extension:
  Xnew = Xnew %*% dmap$X %*% diag(1/dmap$eigenvals)

  return(Xnew)

}
