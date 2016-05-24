rm(list=ls())
library(simone)
set.seed(999)

## Networks generation - the ancestor and its child
p <- 100
ancestor <- rNetwork(p, pi=p, name="ancestor")
child1   <- coNetwork(ancestor, 1, name = "child 1")
child2   <- coNetwork(ancestor, 1, name = "child 2")
## compare the network pairwisely
plot(ancestor, child1)
plot(ancestor, child2)
plot(child1,child2)

## multiple sample setup
n <- c(3*p,3*p)
data  <- rTranscriptData(n,child1,child2)
attach(data)

## Running multitak simone with coopLasso coupling (default)
control <- setOptions(penalty.max=0.1, penalty.min=0.015) 
res     <- simone(X, tasks=tasks, control=control)

## Plotting the results
plot(res, ref.graph=list(child1$Theta,child2$Theta), ask=FALSE)
glist <- getNetwork(res,"BIC")
plot(glist[[1]],glist[[2]])
plot(child1, glist[[1]])
plot(child2, glist[[2]])

detach(data)
