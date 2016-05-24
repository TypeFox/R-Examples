cabe <-
function(dx){
  nf=length(dx)
  nr=sapply(dx,length)
  r1=cumprod(nr)
  repes=r1[length(r1)]/r1
  cabf=mapply(dx,repes,r1[length(r1)],FUN=function(x,y,z) rep(rep(x,each=y),z/(length(x)*y)))
  return(cabf)
}
