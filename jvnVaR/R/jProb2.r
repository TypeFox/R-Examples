jProb2 <-
function(n,k,alpha){
k0 <- floor(alpha * (n+1))
h <- abs(k0-k)
k1 <- max(0, k0-h)
k2 <- min(n, k0+h)
temp <- (n-k0+1)/k0 * alpha/(1-alpha)
p0 <- jBinomial(n,alpha)
p <- p0/ temp
T <- 0
for (i in k0:k2){
temp <- (n-i+1)/i * alpha/(1-alpha)
# temp <- jRatio(n,i,alpha)
p <- p*temp
T <- T + p
}
temp <- (n-k0)/(k0+1) * alpha/(1-alpha)
p <- p0*temp
for (i in k0:k1){
temp <- (n-i)/(i+1) * alpha/(1-alpha)
# temp <- jRatio(n,i+1,alpha)
p <- p/temp
T <- T + p
}
T <- T - p0
return(T/n)
}