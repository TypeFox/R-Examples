## ----echo=FALSE----------------------------------------------------------
  cores <- as.integer(Sys.getenv('SG_RUN'))
  if(is.na(cores)) options(lfe.threads=1)

## ------------------------------------------------------------------------
library(lfe)
set.seed(42)
x <- rnorm(100000)
f1 <- sample(10000,length(x),replace=TRUE)
f2 <- sample(10000,length(x),replace=TRUE)
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)

## ------------------------------------------------------------------------
system.time(est <- felm(y ~ x | f1 + f2))

## ------------------------------------------------------------------------
system.time(alpha <- getfe(est))

## ------------------------------------------------------------------------
f2 <- sample(300,length(x),replace=TRUE)
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x | f1 + f2))
system.time(alpha <- getfe(est))

## ------------------------------------------------------------------------
f2 <- (f1 + sample(10,length(x),replace=TRUE)) %% 300
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x | f1 + f2))
system.time(alpha <- getfe(est))

## ------------------------------------------------------------------------
system.time(est <- felm(y ~ x + factor(f2) | f1))
system.time(alpha <- getfe(est))

## ------------------------------------------------------------------------
f2 <- (f1 + sample(10,length(x),replace=TRUE)^3) %% 300
y <- x + cos(f1) + log(f2+1) + rnorm(length(x), sd=0.5)
system.time(est <- felm(y ~ x | f1 + f2))
system.time(alpha <- getfe(est))
nlevels(est[['cfactor']])

## ----cache=TRUE, tidy=FALSE----------------------------------------------
library(igraph)
mkgraph <- function(f1,f2)
  graph.edgelist(cbind(paste('f1',f1),paste('f2',f2)), directed=FALSE)

appxdiam <- function(g) max(shortest.paths(g,sample(V(g),10),sample(V(g),10)))
f2 <- sample(10000,length(x),replace=TRUE)
appxdiam(mkgraph(f1,f2))
f2 <- (f1 + sample(5,length(x),replace=TRUE)^3) %% 300
appxdiam(mkgraph(f1,f2))
f2 <- (f1 + sample(5,length(x),replace=TRUE)) %% 300
appxdiam(mkgraph(f1,f2))

