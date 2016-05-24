p.WR <- function(N, m, pk){
p <- rep(0,N)
I <- nk(N,m)
N <- dim(I)[1]
for(i in 1:N){
ni <- c(I[i,])
p[i] <- dmultinom(ni, prob=pk)
}
p
}