`rm4pl` <-
function(N=100,S=0,C=0,D=0,s=1/1.702,b=0,c=0,d=1) {
 theta <- data.frame(theta = log(abs(((d-D)-(C+c))/(runif(N)-(C+c)) - 1))/-sqrt(1/(s^2 + S^2)) + b)
 return(theta)
 }

