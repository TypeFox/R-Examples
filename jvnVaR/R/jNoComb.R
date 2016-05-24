jNoComb <-
function(n,k,alpha){
# Number of combination of k elements from a n-set
C <- 1
# n!/k!(n-k)! * alpha^k * (1-alpha)^(n-k)
for (i in 1:k){
#C <- C * (n-i+1)/i
C <- C * (n-i+1)/i*alpha*(1-alpha)
}
C <- C * (1-alpha)^(n-k-k)
return(C)
}
