`ggr4pl` <-
function(n=5, rep=1,theta=0,a=rep(1,n),b=rep(0,n),c=rep(0,n),d=rep(1,n)) {
 result <- NULL
 N      <- rep
 for (j in theta) {
  temp   <- NULL
  for (i in 1:n) temp <- cbind(temp, gr4pl(N=N,theta=j,a[i],b[i],c[i],d[i]))
  result <- rbind(result, temp)
  }
 return(data.frame(result))
 }

