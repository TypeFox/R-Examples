plot.RWBP <-
function(x,...){
V(x$igraph)$color<-ifelse(V(x$igraph)$name<0, 'blue', 'red')
plot(x$igraph, edge.label=round(E(x$igraph)$weight, 3))
#plot(x$igraph)
}
