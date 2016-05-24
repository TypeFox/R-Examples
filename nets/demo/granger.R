
set.seed(1234)

# parameters
N <- 20
T <- 1000
P <- 3

A <- array(0,dim=c(N,N,P))

A[3,7,1] <- 0.2
A[4,1,1] <- 0.2
A[9,1,2] <- 0.2
A[2,5,3] <- 0.2
A[8,3,3] <- 0.2

y      <- matrix(0,T,N)
mu     <- matrix(0,T,N)
eps    <- matrix(rnorm(T*N),T,N)

for( t in (P+1):T ){
  for( l in 1:P ){
    mu[t,] <- mu[t,] + A[,,l] %*% y[t-l,]
  }
  y[t,] <-  mu[t,] + eps[t,]
}

matplot(y,t='l')

# estimate var
lambda  <- 0.01
system.time( mdl <- nets(y,CN=FALSE,p=P,lambda=lambda*T,verbose=TRUE) )

g.adj.hat <- mdl$g.adj

granger.network.hat <- graph.adjacency( g.adj.hat , mode='directed' )

degree <- degree(granger.network)
V( granger.network )$size <- round( (degree/max(degree))*10+2 )

plot( granger.network.hat , edge.arrow.size=0.25 )
