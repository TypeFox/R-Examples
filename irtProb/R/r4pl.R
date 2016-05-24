`r4pl` <-
function(N=100,a=1,b=0,c=0,d=1) {
 theta <- log(abs((d-c)/(runif(N)-c) - 1))/-a + b
 return(theta)
 }
