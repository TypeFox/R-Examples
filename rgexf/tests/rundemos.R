library(rgexf)

demo(edge.list)
demo(gexf)
demo(gexfattributes)
demo(gexfbasic)
demo(gexfbuildfromscratch)
demo(gexfdynamic)
demo(gexfdynamicandatt)
demo(gexffull)
demo(gexftwitter)
demo(gexfigraph)

g <- graph.ring(10)
g <- set.graph.attribute(g, "name", "RING")

# Colors
g <- set.vertex.attribute(g, "color", value=c("red", "green"))

# Weight
g <- set.edge.attribute(g, "weight", value=runif(ecount(g)))

g2 <- igraph.to.gexf(g)
#plot(g)
#plot(g2)
g3 <- gexf.to.igraph(g2)
