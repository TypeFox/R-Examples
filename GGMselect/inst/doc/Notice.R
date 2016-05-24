### R code from vignette source 'Notice.Rnw'

###################################################
### code chunk number 1: exFast
###################################################
 library("GGMselect")
p=30
n=30
# ----------------------------------------
# Random graph generator: use of simulateGraph
# ----------------------------------------
eta=0.11

Gr <- simulateGraph(p,eta)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)
# ----------------------------------------
# Graph selection with family C01:  use of selectFast
# ----------------------------------------
GRest <- selectFast(X, family="C01")



###################################################
### code chunk number 2: Notice.Rnw:787-792
###################################################
# ----------------------------------------
# Plot the result with the help of the package network
# ----------------------------------------
library(network)



###################################################
### code chunk number 3: explotFast
###################################################
gV <- network(Gr$G)
g <- network(GRest$C01$G)
par(mfrow=c(1,2), pty = "s")
a <- plot(gV, usearrows = FALSE)
title(sub="Simulated graph")
plot(g, coord=a, usearrows = FALSE)
title(sub="Graph selected with C01 family")



###################################################
### code chunk number 4: exQE
###################################################
# ----------------------------------------
# Graph selection with family QE:  use of selectQE
# ----------------------------------------
GQE <- selectQE(X)
# ----------------------------------------
# Plot the result
# ----------------------------------------



###################################################
### code chunk number 5: explotQE
###################################################
# CACHER
g <- network(GQE$G)

par(mfrow=c(1,2), pty = "s")
plot(gV,coord=a, usearrows = FALSE)
title(sub="Simulated graph")
plot(g,coord=a, usearrows = FALSE)
title(sub="Graph selected with QE family")



###################################################
### code chunk number 6: exMyFam
###################################################
# ----------------------------------------
# Graph selection with selectMyFam
# ----------------------------------------
# generate a family of candidate graphs with glasso
library("glasso")
MyFamily <- NULL
for (j in 1:3){
  MyFamily[[j]] <- abs(sign(glasso(cov(X),rho=j/5)$wi))
  diag(MyFamily[[j]]) <- 0
}
# select a graph within MyFamily
GMF <- selectMyFam(X,MyFamily)
# ----------------------------------------
# Plot the result
# ----------------------------------------




###################################################
### code chunk number 7: explotMyFam
###################################################
# CACHER
# plot the result
g <- network(GMF$G)

par(mfrow=c(1,2), pty = "s")
plot(gV,coord=a, usearrows = FALSE)
title(sub="Simulated graph")
plot(g,coord=a, usearrows = FALSE)
title(sub="Graph selected with MyFam")
# FIN CACHER


