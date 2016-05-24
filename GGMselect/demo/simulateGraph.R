# simulate a graph
p=30
eta=0.13
Gr <- simulateGraph(p,eta)

# plot the graph
library(network)
gV <- network(Gr$G)
plot(gV,jitter=TRUE, usearrows = FALSE, label=1:p,displaylabels=TRUE)

