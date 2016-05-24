# +++++++++++++++++++++++++++++++++++++++++++++++++
p=30
n=30
# simulate graph
eta=0.11
Gr <- simulateGraph(p,eta)
# simulate data
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)
# estimate graph
GQE <- selectQE(X)

# plot the result
library(network)
old.par <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
gV <- network(Gr$G)
a<-plot(gV, usearrows = FALSE)
gQE <- network(GQE$G)
plot(gQE, usearrows = FALSE, coord=a)
par(old.par)


