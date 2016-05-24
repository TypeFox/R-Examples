library(PEIP)
set.seed(11)
####  perfect data with no noise
n <- 5
A <- matrix(runif(n*n),nrow=n)
B <- runif(n)
###  get right-hand-side (data)
trhs = as.vector( A %*% B  )
Lout = rlsqr(G=A,u=trhs ,wantse=0,damp=0,
    atol=0,btol=0,conlim=0,itnlim=100,nout=0)

Lout$x
B

sqrt(  sum( (Lout$x - B)^2 ) /n )


