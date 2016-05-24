#' calculate the number of tips per origin (NoTO) metric
#'
#' This function calculates the NoTO metric. To conduct significance test on NoTO, use 'noto' as input of 'func' in 'treestat'
#' @param phy an object of class 'phylo'.
#' @param state a vector of '0' and '1' for trait state of each tip, in the same order as the tip labels.

noto <- function (state,phy) {
	nspecies<-length(state)
	aa<-order(phy$edge[,1],decreasing=T)
	phy$edge<-phy$edge[aa,]
	tippos <- which(phy$edge[,2] %in% which(state==1))
	anc<-phy$edge[seq(from=1,by=2,length.out=length(phy$edge[,1])/2),1]
	des<-matrix(phy$edge[,2],ncol=2,byrow=T)
	I<-rep(NA,phy$Nnode)
	I<-c(state,I)
	I<-cbind(I,I)
	med<-function(x) {
		i<-length(x)/2
		sort(x)[c(i,i+1L)]
	}
	for (i in 1:phy$Nnode) {
		I[anc[i],]<-med(I[des[i,],])
		}
	nodestate<-matrix(NA,(phy$Nnode*2+1),2)
	Iw<-as.vector(I[des[phy$Nnode,],])
	nodestate[(phy$Nnode+2),]<-range(med(c(0,0,Iw)))
	for (i in (phy$Nnode-1):1) {
		j<-anc[i]
		Iw<-as.vector(I[des[i,],])
		k<-which(phy$edge[,2]==j)
		tmp<-nodestate[phy$edge[k,1],]
		nodestate[j,1]<-min(med(c(tmp[1],tmp[1],Iw)))
		nodestate[j,2]<-max(med(c(tmp[2],tmp[2],Iw)))
	}
	un<-which(nodestate[,1]!=nodestate[,2])
	nodestate[un,1]<-as.numeric(runif(length(un))<=0.5)
	nodestate<-nodestate[,1]
	origin<-numeric()
	for(k in tippos) {
		parent<-which(phy$edge[,2]==phy$edge[k,1])
		if((length(parent)==0)||(nodestate[phy$edge[parent,2]]!=1)) {
			origin<-c(origin,k)
		}else {
			while(nodestate[phy$edge[parent,2]]==1) {
				parent2<-which(phy$edge[,2]==phy$edge[parent,1])
				if((length(parent2)==0)||(nodestate[phy$edge[parent2,2]]!=1)) {
					break
				}else {
					parent<-parent2
				}
			}
		origin<-c(origin,parent)
		}
	}
	sum(state)/length(unique(phy$edge[origin,1]))
}
