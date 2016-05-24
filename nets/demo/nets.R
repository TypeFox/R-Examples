
library(MASS)
library(nets)

set.seed(1234)

N <- 20
P <- 2
L <- 40
T <- 1000

A <- array(0,dim=c(N,N,P))
C <- matrix(0,N,N)

A[,,1]   <- 0.7 * diag(N) 
A[,,2]   <- 0.2 * diag(N) 
A[1,2,1] <- 0.2
A[4,3,2] <- 0.2

C      <- diag(N)
C[1,1] <- 2
C[4,2] <- -0.2
C[2,4] <- -0.2
C[1,3] <- -0.1
C[1,3] <- -0.1

Sig    <- solve(C)

y      <- matrix(0,T,N)
eps    <- mvrnorm(T,rep(0,N),Sig)

for( t in (P+1):T ){
  for( l in 1:P ){
    y[t,] <- y[t,] + A[,,l] %*% y[t-l,]
  }
  y[t,] <-  y[t,] +  eps[t,]
}

#
lambda <- c(0.1,0.01)
system.time( mdl <- nets(y,p=P,lambda=lambda*T,verbose=TRUE)  )

g.adj.hat <- mdl$g.adj

granger.network.hat <- graph.adjacency( g.adj.hat , mode='directed' )

degree <- degree(granger.network)
V( granger.network )$size <- round( (degree/max(degree))*10+2 )

plot( granger.network.hat , edge.arrow.size=0.25 )

c.adj.hat <- mdl$c.adj
contemporaneous.network <- graph.adjacency( c.adj.hat , mode='undirected' )

degree <- degree(contemporaneous.network)
V( contemporaneous.network )$size <- round( (degree/max(degree))*10+2 )

plot( contemporaneous.network )
