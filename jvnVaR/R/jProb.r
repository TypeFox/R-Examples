jProb <-
function(n,k1,k2,alpha){
mu <- n*alpha
sigma <- sqrt(n*alpha*(1-alpha))
T <- jSimson(1000,(k1-mu)/sigma,(k2-mu)/sigma,jNormPdf)
return(T)
}
