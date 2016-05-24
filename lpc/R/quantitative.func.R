quantitative.func  <- function(x,y,s0=0){
  # regression of x on y
  my=mean(y)
  yy <- y-my
  temp <- x%*%yy
  mx=rowMeans(x)
  syy= sum(yy^2)
  scor <- temp/syy
  b0hat <- mx-scor*my
  ym=matrix(y,nrow=nrow(x),ncol=ncol(x),byrow=T)
  xhat <- matrix(b0hat,nrow=nrow(x),ncol=ncol(x))+ym*matrix(scor,nrow=nrow(x),ncol=ncol(x))
  sigma <- sqrt(rowSums((x-xhat)^2)/(ncol(xhat)-2))
  sd <- sigma/sqrt(syy)
  tt <- scor/(sd+s0)
  return(list(tt=tt, numer=scor, sd=sd))
}
