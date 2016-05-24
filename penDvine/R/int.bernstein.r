int.bernstein <- function(penden.env,Y) {
  int.base <- matrix(0,length(Y),(get("K",penden.env)+1))
  index.k <- matrix(1:length(Y))
  for(j in 1:(get("K",penden.env)+1)) int.base[,j] <- spline(get(x="x.help",penden.env),y=get("int.base",penden.env)[,j],xout=Y)$y
  return(int.base)
}
