plottanet <-
function(timeaggregatednetwork,layout=layout.circle,vertex.label=V(timeaggregatednetwork)$name,vertex.size=0,vertex.label.cex=0.5,edge.arrow.size=0.5,edge.width=E(timeaggregatednetwork)$Count/5)
{
	plot(timeaggregatednetwork,layout=layout,vertex.label=vertex.label,vertex.size=vertex.size,vertex.label.cex=vertex.label.cex,edge.arrow.size=edge.arrow.size,edge.width=edge.width,vertex.label.family="Helvetica")	
}

