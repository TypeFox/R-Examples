ews.trans <-
function(x,scales=NULL){
  N=length(x)
  if (length(scales)==0) J=floor(log(N,2))
  else {J=length(scales)}
  res=wavMODWT(x=x, wavelet="haar")
  res=as.matrix(res)
  out=matrix(res[1:(J*N)],N,J)
  for (i in 1:(J)) out[,i]=c(out[(2^i):N,i],out[1:(2^i-1),i])
  return(out^2)
}
