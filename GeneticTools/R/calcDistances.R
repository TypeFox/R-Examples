# Experimental parts on time series clustering and the corresponding distances measures

calcDistances <- function(X,method="C",nodes=20){
  res <- c()
  if(method=="C") res <- calcDistances.C(X)
  if(method=="splines") res <- calcDistances.Splines(X,nodes=nodes)
  res
} 

calcDistances.Splines <- function(X,nodes){
  NC <- ncol(X)
  NR <- nrow(X) 
  for(i in 1:ncol())
  {
    temp <- smooth.spline(1:NC, X[i,], df=NC, all.knots=TRUE)
    deriv0 <- predict(temp,seq(1,NC,0.1),deriv=0)$y
    deriv1 <- predict(temp,seq(1,NC,0.1),deriv=1)$y
  }
}