
#' @keywords internal

auxmat <-
function(m,vettcsi,vettphi,a,b,c,d,e){
  elemat<-matrix(NA,nrow=m,ncol=length(vettcsi))
  for(k in 1:m){
    elemat[k,]<-e*((k-1)^d)/((a+b*vettcsi+vettphi*(k-1))^c)
  }
  return(elemat)
}


