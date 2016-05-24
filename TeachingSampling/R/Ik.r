Ik <- function(N,n){
Q <- Support(N,n,ID=FALSE)
I <- matrix(0,choose(N,n),N)
for(i in 1:n){
for(j in 1:choose(N,n)){
for(k in 1:N){
if (Q[j,i]==k)
I[j,k] <- 1
}
}
}
I
}