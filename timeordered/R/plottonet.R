plottonet <-
function(g,path=NULL,edgecolor="gray",edgehighlightcolor="red",vertex.size=0.01,edge.arrow.size=0.1,edge.width=0.2,vertex.color=NA,vertex.label.cex=0.1,vertex.frame.color=NA,vertex.label.color="black")
{
	E(g)$color=edgecolor
	if (length(path) > 0)
	{
		E(g, path=path)$color=edgehighlightcolor
	}
	names <- V(g)$Name
	names[V(g)$Time > min(V(g)$Time)] <- NA

	plot(g, layout=cbind(as.numeric(factor(V(g)$Name)), V(g)$Time), vertex.label=names,vertex.size=vertex.size,edge.arrow.size=edge.arrow.size,edge.width=edge.width, 	vertex.color=vertex.color,vertex.label.cex=vertex.label.cex,vertex.frame.color=vertex.frame.color,vertex.label.color=vertex.label.color,vertex.label.family="Helvetica")
}

