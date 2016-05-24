"sackin" <-
function(tree, norm=NULL){ 
if (identical(tree,NULL)) {
		stop("invalid tree","\n")
		
}
	clades=smaller.clade.spectrum(tree)
	INS=sum(clades[,1])
	
	
	if (identical(norm,NULL)==TRUE) { return(INS) }
	if (norm=="pda") {
		leaf_nb<-nrow(tree$merge) + 1
		res <- INS/((nrow(tree$merge) + 1)^(3/2))
		return(res)
	}
	if (norm=="yule") {
		leaf_nb<-nrow(tree$merge) + 1
		EINS <- 2*leaf_nb*sum(1/2:leaf_nb)
		res <- (INS-EINS)/(leaf_nb)
		return(res)
	}
	stop("Incorrect argument for 'norm'")
	 
}

