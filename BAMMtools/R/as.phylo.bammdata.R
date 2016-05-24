as.phylo.bammdata <- function(x, ...) {
	
	if (!'bammdata' %in% class(x)) {
		stop("Object ephy must be of class bammdata\n");
	}		
	
	newphylo <- list();
	newphylo$edge <- x$edge;
	newphylo$Nnode <- x$Nnode;
	newphylo$tip.label <- x$tip.label;
	newphylo$edge.length <- x$edge.length;
	class(newphylo) <- 'phylo';
	attributes(newphylo)$order = attributes(x)$order;
	if (attributes(newphylo)$order != "cladewise") {
		newphylo <- reorder(newphylo);
	}
	return(newphylo);
}
