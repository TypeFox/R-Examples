my.bisec <-
function(f,int.l,int.u,eps=1e-6){
  while(max((int.l-int.u)^2)>eps){
    int.m <- (int.l+int.u)/2
    fl <- f(int.l)
    fu <- f(int.u)
    fm <- f(int.m)
    signs <- sign(fl*fu)
	i0 <- which(signs>=0)
    if (length(i0)>0){
      int.l[i0] <- ifelse(abs(fl[i0])>=abs(fu[i0]),int.u[i0],int.l[i0])
      int.u[i0] <- ifelse(abs(fl[i0])>=abs(fu[i0]),int.u[i0],int.l[i0])
    }
	i2 <- which(signs==-1)
    if (length(i2)>0){
    int.l[i2] <- ifelse(sign(fl[i2])==sign(fm[i2]),int.m[i2],int.l[i2])
    int.u[i2] <- ifelse(sign(fl[i2])==sign(fm[i2]),int.u[i2],int.m[i2])
	}
  }
  return(int.l)
}
