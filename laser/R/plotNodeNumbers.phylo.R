`plotNodeNumbers.phylo` <-
function(phy)
{
	phy$node.label <- (length(phy$tip.label)+1):max(phy$edge);
	plot.phylo(phy, show.node.label=TRUE);	
}

