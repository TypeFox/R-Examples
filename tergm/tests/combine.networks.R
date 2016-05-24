# This tests the internal function combine.networks.

n <- 10
T <- 10

library(tergm)
yl <- replicate(T,
                {
                  y <- network.initialize(n,dir=FALSE)
                  y <- simulate(y~edges, coef=-1)
                  y %n% "x" <- matrix(runif(n*n),n,n)
                  y %v% "v" <- runif(n)
                  y %e% "e" <- runif(network.edgecount(y))
                  y
                },
                simplify=FALSE)

yc <- tergm:::combine.networks(yl)
ym <- as.matrix(yc)

for(t in seq_len(T)){
  J <- I <- (t-1)*n + seq_len(n)
  stopifnot(all(as.matrix(yc)[I,J]==as.matrix(yl[[t]]))) # Check ties.
  stopifnot(all((yc %n% "x")[I,J]==(yl[[t]] %n% "x"))) # Check dyadic attributes.
  stopifnot(all((yc %v% "v")[I]==(yl[[t]] %v% "v"))) # Check vertex attributes.
  stopifnot(all(as.matrix(yc,attr="e")[I,J]==as.matrix(yl[[t]],attr="e"))) # Check edge attributes.

  ym[I,J] <- 0
}
stopifnot(all(ym==0))


m <- 7

yl <- replicate(T,
                {
                  y <- network.initialize(n,dir=FALSE, bipartite=m)
                  y <- simulate(y~edges, coef=-0.5)
                  y %n% "x" <- matrix(runif(m*(n-m)),m,n-m)
                  y %v% "v" <- runif(n)
                  y %e% "e" <- runif(network.edgecount(y))
                  y
                },
                simplify=FALSE)

yc <- tergm:::combine.networks(yl)
stopifnot(identical(yc%n%"bipartite",T*m))

ym <- as.matrix(yc)

for(t in seq_len(T)){
  I <- (t-1)*m + seq_len(m)
  J <- (t-1)*(n-m) + seq_len(n-m)
  stopifnot(all(as.matrix(yc)[I,J]==as.matrix(yl[[t]]))) # Check ties.
  stopifnot(all((yc %n% "x")[I,J]==(yl[[t]] %n% "x"))) # Check dyadic attributes.
  stopifnot(all((yc %v% "v")[I]==(yl[[t]] %v% "v")[1:m])) # Check vertex attributes for egos.
  stopifnot(all((yc %v% "v")[T*m+J]==(yl[[t]] %v% "v")[(m+1):n])) # Check vertex attributes for alters.
  stopifnot(all(as.matrix(yc,attr="e")[I,J]==as.matrix(yl[[t]],attr="e"))) # Check edge attributes.

  ym[I,J] <- 0
}
stopifnot(all(ym==0))


