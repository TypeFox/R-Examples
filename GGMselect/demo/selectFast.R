# +++++++++++++++++++++++++++++++++++++++++++++++++
p=30
n=30
# simulate graph
eta=0.11
Gr <- simulateGraph(p,eta)
X <- rmvnorm(n, mean=rep(0,p), sigma=Gr$C)
# estimate graph
GRest <- selectFast(X, family="C01")

# plot result
library(network)
old.par <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
gV <- network(Gr$G)
a <- plot(gV, usearrows = FALSE)
g <- network(GRest$C01$G)
plot(g,coord=a, usearrows = FALSE)
par(old.par)


