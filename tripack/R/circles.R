circles <- function(x,y,r,...){
  n <- length(x)
  if(length(y)!=n || length(r)!=n)
    stop("arguments should be of same length!")
  phi <- seq(0,2*pi,length=360)
  for(i in 1:n){
    lines(c(x[i]+r[i]*cos(phi),x[i]+r[i]),c(y[i]+r[i]*sin(phi),y[i]),type="l",...)
  }
}
