#' calculate the tip age rank sum (TARS) metric
#'
#' This function calculates the TARS metric. To conduct significance test on TARS, use 'tars' as input of 'func' in 'treestat'
#' @param phy an object of class 'phylo'.
#' @param state a vector of '0' and '1' for trait state of each tip, in the same order as the tip labels.

tars <- function (state,phy) {
	nspecies<-length(state)
	aa<-order(phy$edge[,1],decreasing=T)
	phy$edge<-phy$edge[aa,]
	phy$edge.length<-phy$edge.length[aa]
	BL<-phy$edge.length[order(phy$edge[,2])][1:nspecies]
	tiprank<-rank(BL)
	sumrank<-sum(tiprank[state==1])
}
