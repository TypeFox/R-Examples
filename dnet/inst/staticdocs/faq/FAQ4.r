# <span style="font-weight:bold; color:#F87217; text-decoration:underline">The dnet package supports two ways to visualise the identified networks/graphs: 1) the network itself as a single display, 2) the same network but with multiple colorings/displays according to samples.</span>

# To demonstrate the visuals supported, we use a random graph generated according to the ER model, and only keep the largest component.
g <- erdos.renyi.game(100, 1/100)
g <- dNetInduce(g, nodes_query=V(g))

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">To display the network itself, the key setting is the layout.</span> Since the network is represented as an object of class 'igraph', the layouts supported in the 'igraph' package can be found in the <a href="http://igraph.org/r/doc/layout.html">layout</a> page. Below is an example using a force-based algorithm proposed by Fruchterman and Reingold:
visNet(g, layout=layout.fruchterman.reingold)
# Additionally, the dnet package also uses two other layouts/diagrams: 1) arc diagram in one-dimensional layout, and 2) circle diagram. 

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">To better display the network, it is advisable to incorporate community information.</span> In doing so, we first identify communities based on a spin-glass model and simulated annealing.
com <- igraph::spinglass.community(g, spins=25)
vgroups <- com$membership
## color nodes: according to communities
mcolors <- visColormap(colormap="rainbow")(length(com))
vcolors <- mcolors[vgroups]
## size nodes: according to degrees
vdegrees <- igraph::degree(g)
## sort nodes: first by communities and then degrees
df <- data.frame(ind=1:vcount(g), vgroups, vdegrees)
ordering <- df[order(vgroups,vdegrees),]$ind

# Now, make comparsions between different layouts
## using 1-dimensional arc diagram
visNetArc(g, ordering=ordering, vertex.label.color=vcolors, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.size=log(vdegrees)+0.1)
## using circle diagram (drawn into a single circle)
visNetCircle(g, colormap="rainbow", com=com, ordering=ordering)
## using circle diagram (drawn into multlpe circles: one circle per community)
visNetCircle(g, colormap="rainbow", com=com, circles="multiple", ordering=ordering)
## using a force-based algorithm proposed by Fruchterman and Reingold
visNet(g, colormap="rainbow", layout=layout.fruchterman.reingold, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere")
## when using force-based layout, it is also useful to highlight the communities in the background, and have edges being colored differently according whether an edge lies within a community or between communities. 
mark.groups <- igraph::communities(com)
mark.col <- visColoralpha(mcolors, alpha=0.2)
mark.border <- visColoralpha(mcolors, alpha=0.2)
edge.color <- c("grey", "black")[igraph::crossing(com,g)+1]
visNet(g, colormap="rainbow", glayout=layout.fruchterman.reingold, vertex.color=vcolors, vertex.frame.color=vcolors, vertex.shape="sphere", mark.groups=mark.groups, mark.col=mark.col, mark.border=mark.border, mark.shape=1, mark.expand=10, edge.color=edge.color)

# Now, let us look at <span style="font-weight:bold; color:#F87217; text-decoration:underline"> visualising the same network but with the multiple colorings/displays according to samples.</span>

# Assume we have 10 samples, each containing numeric information about nodes in the graph:
nnodes <- vcount(g)
nsamples <- 10
data <- matrix(runif(nnodes*nsamples), nrow=nnodes, ncol=nsamples)
rownames(data) <- V(g)$name

# <span style="font-weight:bold; color:#F87217; text-decoration:underline">In the dnet package, these sample-specific network visuals can be simply as they are provided, or samples being self-organised onto 2D landscape.</span>

## simply as they are provided
visNetMul(g, colormap="rainbow", data=data, glayout=layout.fruchterman.reingold)

## being self-organised onto 2D sample landscape (ie a sheet-shape rectangle grid)
sReorder <- dNetReorder(g, data, feature="node", node.normalise="none")
visNetReorder(g, colormap="rainbow", data=data, sReorder)

## By default, 2D sample landscape is built based on node features without considering node degrees. To take into account the connectivity in the network, the information used can be on edges which are transformed from information on nodes: input data and degree. In doing so, <span style="font-weight:bold; color:#F87217; text-decoration:underline">the transformed matrix of network edges × samples are used for self-organising samples onto 2D landscape (implemented in the 'supraHex' package).</span>
sReorder <- dNetReorder(g, data, feature="edge", node.normalise="degree")
visNetReorder(g, colormap="rainbow", data=data, sReorder)
