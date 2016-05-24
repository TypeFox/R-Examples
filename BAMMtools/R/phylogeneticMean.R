phylogeneticMean <- function(traits, phy, lambda=1){
	
	if (!is.null(names(traits)))
		traits <- traits[phy$tip.label];
	
	if (class(phy) == 'phylo'){
		vmat <- vcv.phylo(phy);			
	}else{
		vmat <- phy;
		
	}
	
	
	dd <- diag(vmat);
	vmat <- vmat * lambda;
	diag(vmat) <- dd;
	
	onev <- matrix(rep(1, length(traits)), nrow=length(traits), ncol=1);

	anc <- as.vector(solve(t(onev)%*% solve(vmat) %*% onev) %*% (t(onev)%*% solve(vmat) %*% traits));
	beta <- as.vector((t(traits-anc) %*% solve(vmat) %*% (traits-anc))/length(traits));
	
	return(list(anc=as.vector(anc), beta=as.vector(beta)));

}
