rm(list=ls())
library(simone)
## Data set and graph generation: a network with no latent
## clustering, number of edges = number of nodes
set.seed(777)
p    <- 50
n    <- 3*p # quite favorable settings :)
g    <- rNetwork(p, pi=p, name="Theorical Network")
data <- rTranscriptData(n,g)
attach(data)

## Running simone
res <- simone(X, type="steady-state")

## Plotting the results
plot(res, ref.graph=g$A, ask=FALSE)

## Compare the network (assuming the number of edges is known)
plot(g,getNetwork(res, p))

detach(data)
