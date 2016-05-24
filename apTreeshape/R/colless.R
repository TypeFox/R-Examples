"colless" <-
function(tree, norm=NULL){

if (identical(tree,NULL)) {
		stop("invalid tree","\n")
}	

norm.yule <- function(ICN, tree) {
	leaf_nb <- nrow(tree$merge) + 1
# 	if (leaf_nb<30){
# 		EICN <- lsmall_mean_coll_yule[leaf_nb] 
# 	} 
# 	else {
#gamma = 0.57721566
		EICN <- leaf_nb*log(leaf_nb) + (0.57721566 -1 - log(2))*leaf_nb
# 	}
	IC <- (ICN - EICN)/(leaf_nb)
	IC
}

norm.pda <- function(ICN, tree) {
	leaf_nb <- nrow(tree$merge) + 1
	IC <- ICN/((leaf_nb)^(3/2))
	IC
}

	
	tmp=smaller.clade.spectrum(tree)
	res=sum(abs(tmp[,1]-2*tmp[,2]))
	
	
	if (identical(norm, NULL)==TRUE) { return(res) }
	if (norm=="pda") { return(norm.pda(res, tree)) }
	if (norm=="yule") { return(norm.yule(res, tree)) }
	
	stop("Incorrect argument for 'norm'")
	
}

