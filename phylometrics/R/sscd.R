#' calculate the sum of sister clade differences (SSCD) metric
#'
#' This function calculates the SSCD metric. To conduct significance test on SSCD, use 'sscd' as input of 'func' in 'treestat'
#' @param phy an object of class 'phylo'.
#' @param state a vector of '0' and '1' for trait state of each tip, in the same order as the tip labels.

sscd <- function(state,phy) {
	aa<-order(phy$edge[,1],decreasing=T)
	phy$edge<-phy$edge[aa,]
	anc<-phy$edge[seq(from=1,by=2,length.out=length(phy$edge[,1])/2),1]
	des<-matrix(phy$edge[,2],ncol=2,byrow=T)
	I<-rep(NA,phy$Nnode)
	I<-c(state,I)
	for (i in 1:phy$Nnode) {
		I[anc[i]]<-mean(I[des[i,]])
	}
	out<-I[des]
	out<-abs(out[1:phy$Nnode]-out[(phy$Nnode+1):(2*phy$Nnode)])
	sum(out)
}
