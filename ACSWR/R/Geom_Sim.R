Geom_Sim <-
function(p,n){
  q <- 1-p
  x <- numeric(n)
  for(i in 1:n){
    temp <- runif(1)
    temp <- 1-temp
    j <- 0
    while(((temp>q^j) & (temp <= q^{j-1}))==FALSE)j <- j+1
    x[i] <- j
  }
  return(x)
}
