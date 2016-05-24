int.bernstein.help <- function(penden.env) {
  len.k <- get("K",penden.env)
  x.help <- seq(0,1,length=501)
  int.base <- matrix(0,length(x.help),(len.k+1))
  for(j in 0:len.k) int.base[,(j+1)] <- apply(matrix(x.help,ncol=1),1,help2,v=j,n=len.k)
  assign("int.base",int.base,penden.env)
  assign("x.help",x.help,penden.env)
}
