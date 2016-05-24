#' calculate the Fritz & Purvis D (FPD) metric
#'
#' This function calculates the FPD metric. To conduct significance test on FPD, use 'fpd' as input of 'func' in 'treestat'
#' @param phy an object of class 'phylo'.
#' @param state a vector of '0' and '1' for trait state of each tip, in the same order as the tip labels.
fpd <- function (state,phy) {
	m<-sum(state)
    nspecies<-length(state)
	aa<-order(phy$edge[,1],decreasing=T)
	phy$edge<-phy$edge[aa,]
	phy$edge.length<-phy$edge.length[aa]
	anc<-phy$edge[seq(from=1,by=2,length.out=length(phy$edge[,1])/2),1]
	des<-matrix(phy$edge[,2],ncol=2,byrow=T)
	DES<-matrix(0,phy$Nnode,nspecies)
	for (i in 1:phy$Nnode) {
		tmp<-des[i,]
		offs<-which(anc %in% des[i,])
		while (length(offs)>0) {
			tmp<-c(tmp,des[offs,])
			offs<-which(anc %in% des[offs,])
		}
		tmp<-tmp[tmp<=nspecies]
		DES[i,tmp]<-1
	}
	DES<-rbind(diag(nspecies),DES)
	BL<-phy$edge.length[c(order(phy$edge[,2])[1:nspecies],order(phy$edge[,2],decreasing=T)[-c((nspecies-1):(2*nspecies-2))])]
	BL<-c(BL,0)
	V<-crossprod(DES,DES*BL)
	bmstate<-mvtnorm::rmvnorm(n=1000,sigma=V)
	bmthred<-apply(bmstate,1,quantile,m/nspecies)
	bmstate<-sweep(bmstate,1,bmthred,'<')
	bmstate<-matrix(as.numeric(bmstate),1000,nspecies)
	rdstate<-replicate(1000,sample(state,size=nspecies,replace=F))
	rdstate<-t(rdstate)
	calD <- function (state,anc,des,phy) {
		I<-rep(NA,phy$Nnode)
		I<-c(state,I)
		for (i in 1:phy$Nnode) {
			I[anc[i]]<-mean(I[des[i,]])
		}
		out<-I[des]
		out<-abs(out[1:phy$Nnode]-out[(phy$Nnode+1):(2*phy$Nnode)])
		sum(out)
	}
	obsD<-calD(state,anc,des,phy)
	bmD<-apply(bmstate,1,calD,anc,des,phy)
	rdD<-apply(rdstate,1,calD,anc,des,phy)
	(obsD-mean(bmD))/(mean(rdD)-mean(bmD))
}
