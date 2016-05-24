"is.binary.phylo" <-
function(phy) {
	
	if (class(phy)!='phylo') {
		stop("invalid arguments")
	}
		
	x<-0
	tmp<-rep(x, nrow(phy$edge)+1)
	for (node in 1:(nrow(phy$edge))) {
		tmp[phy$edge[node,1]]<-tmp[phy$edge[node,1]]+1
	}

	res1<-seq(1)
	res2<-seq(1)
	for (node in 1:(nrow(phy$edge)+1)) {
		if (!(tmp[node]==2 | tmp[node]==0)) {
			res1<-c(res1, node)
			res2<-c(res2, tmp[node])
		}
	}
	
	if (length(res1)==1) {
		TRUE
	}
	else {
		res1<-res1[-1]
		res2<-res2[-1]
#		cat("Arbre non binaire :", length(res1), "noeud(s) non binaire(s)\n")
		res=matrix(0,nrow=length(res1),ncol=2)
		res[,1]<-res1
		res[,2]<-res2
		res
	}
}

