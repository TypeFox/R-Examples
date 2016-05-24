rm(list=ls())
library(simone)
set.seed(777)

##-------------------------------
## SIMULATIONS : Data set and graph generation
p     <- 50
n     <- 150
alpha <- c(0.1,0.9)
pi    <- t(matrix(c(.05,.3,0,0),2,2))
g     <- rNetwork(p, pi = pi, alpha = alpha, directed=TRUE)
data  <- rTranscriptData(n,g)

control <- setOptions(edges.max=150)

attach(data)
res.no <- simone(X, type="time-course", control=control)
g.no <- getNetwork(res.no, 80)
plot(g,g.no)

res.cl <- simone(X, type="time-course", clustering=TRUE, control=control)
g.cl <- getNetwork(res.cl, 80)
plot(g,g.cl)

g.no$name <- "without clustering"
g.cl$name <- "with clustering"
plot(g.no,g.cl)

detach(data)
