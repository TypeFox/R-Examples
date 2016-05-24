
library(MASS)

set.seed(1234)

N <- 50
T <- 500

C    <- diag(N)*0.5 + as.matrix( graph.laplacian(erdos.renyi.game(N,N,type='gnm'),normalized=TRUE) )
Sig  <- solve(C)
y    <- mvrnorm(T,rep(0,N),Sig)

lambda  <- 0.03
system.time( mdl <- nets(y,GN=FALSE,lambda=lambda*T,verbose=TRUE) )

c.adj.hat <- mdl$c.adj
contemporaneous.network <- graph.adjacency( c.adj.hat , mode='undirected' )

degree <- degree(contemporaneous.network)
V( contemporaneous.network )$size <- round( (degree/max(degree))*10+2 )

plot( contemporaneous.network )
