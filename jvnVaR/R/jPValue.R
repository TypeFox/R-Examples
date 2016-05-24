jPValue <-
function(n,k,alpha,test_significant){
# n: number of observation
# k: number of outstanding
# alpha: null-hypothesis probability of outstanding
if (n>1000000){
k0 <- floor(alpha * (n+1))
h <- abs(k0-k)
k1 <- max(0, k0-h)
k2 <- min(n, k0+h)
p <- jProb(n,k1,k2,alpha)
}
if (n<=1000000){
p <- jProb2(n,k,alpha)
}
p_value <- 1-p
rslt <- p_value >= test_significant
return(c(p_value,test_significant,rslt))
}
