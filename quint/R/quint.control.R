quint.control <-
function (crit="es",maxl=10,a1=NULL, a2=NULL, w=NULL,Bootstrap=TRUE,B=25,dmin=0.30) 
{ #crit="es" (effect size criterion) or "dm" (difference in means criterion)
  #maxl: maximum total number of leaves (terminal nodes) of the final tree: Lmax
  
  if(!is.null(a1)){
    parvec=round(c(a1,a2))}
  else{parvec=NULL}
  list(crit=crit,maxl=maxl,parvec=parvec,w=w,Boot=Bootstrap,B=B,dmin=dmin)
}
