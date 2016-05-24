WordCountAna<-function(base,sep.word=NULL,ncp=Inf,nb.panel=3,nb.simul=500,proba=0.05,graph=TRUE,axes=c(1,2)){
	within.inertia<-function(coord,weight) sum(t(t(coord)-apply(coord*weight,2,sum)/sum(weight))^2*weight)
	base<-as.data.frame(base)
	if (!inherits(base, "data.frame")) stop("base should be a data.frame")
	groups<-numeric()
	ind.words<-numeric()
	data <- matrix(0,nrow(base),0)
	for (i in 1:ncol(base)){
		ind.table<-textual(tab=cbind.data.frame(rownames(base),base[,i]),num.text=2,contingence.by=1,maj.in.min = TRUE,sep.word=sep.word)$cont.table
		ind.words<-c(ind.words,colnames(ind.table))
		colnames(ind.table)<-paste(colnames(ind.table),i,sep="")
		groups[i]<-ncol(ind.table)
		data<-cbind.data.frame(data,as.data.frame(ind.table))
	}
	dist.words<-levels(as.factor(ind.words))
	size.dist.words<-sapply(dist.words,function(dist.word,ind.words) sum(data[,which(ind.words==dist.word)]),ind.words)
	length.dist.words<-sapply(dist.words,function(dist.word,ind.words) length(which(ind.words==dist.word)),ind.words)
	info.dist.words<-as.data.frame(matrix(ncol=2,nrow=length(dist.words),dimnames=list(dist.words,c("Nb times","Nb panellists"))))
	info.dist.words[,1]<-size.dist.words
	info.dist.words[,2]<-length.dist.words
	info.dist.words<-info.dist.words[order(info.dist.words[,2],decreasing = TRUE),]
	res.mfact<-MFA(data,type=rep("f",ncol(base)),group=groups,ncp=ncp,graph=FALSE,name.group=1:ncol(base))
	sel.words<-names(which(sapply(dist.words,function(dist.word,ind.words) length(which(ind.words==dist.word)),ind.words)>=nb.panel))
	within<-sapply(sapply(sel.words,function(sel.word,ind.words) which(ind.words==sel.word),ind.words),function(pos) within.inertia(res.mfact$freq$coord[pos,],res.mfact$global.pca$call$col.w[pos]))
	within.dist<-list()
	size.sel.words<-sapply(sel.words,function(sel.word,ind.words) sum(data[,which(ind.words==sel.word)]),ind.words)
	length.sel.words<-sapply(sel.words,function(sel.word,ind.words) length(which(ind.words==sel.word)),ind.words)
	levels.length<-levels(as.factor(length.sel.words))
	for(i in 1:length(levels.length)){
		Rsamples<-replicate(nb.simul,sample(1:nrow(res.mfact$freq$coord),levels.length[i],replace=TRUE))
		within.dist[[i]]<-apply(Rsamples,2,function(pos) within.inertia(res.mfact$freq$coord[pos,],res.mfact$global.pca$call$col.w[pos]))
	}
	names(within.dist)<-levels.length
	pvalues<-numeric()
	for (i in 1:length(sel.words)) pvalues[i]<-length(which(within.dist[[which(names(within.dist)==length.sel.words[i])]]<within[i]))/nb.simul
	consensus<-as.data.frame(matrix(ncol=3,nrow=length(sel.words),dimnames=list(sel.words,c("Nb times","Nb panellists","Pvalue"))))
	consensus[,1]<-size.sel.words
	consensus[,2]<-length.sel.words
	consensus[,3]<-pvalues
	consensual.words<-which(consensus[,3]<proba)
	consensus<-consensus[order(consensus[,3]),]
	centroids<-t(sapply(sapply(dist.words,function(dist.word,ind.words) which(ind.words==dist.word),ind.words),function(pos) apply(sweep(as.data.frame(res.mfact$freq$coord)[pos,],1,res.mfact$global.pca$call$col.w[pos],"*"),2,sum)/sum(res.mfact$global.pca$call$col.w[pos])))
	resultats<-list(mfact=res.mfact,dist.words=info.dist.words,centroids=centroids,cons=consensus,cons.words=sel.words[consensual.words])
	class(resultats) <- c("WordCountAna", "list")
	if (graph){
		plot.WordCountAna(resultats,choix="prod")
		plot.WordCountAna(resultats,choix="panel")
		plot.WordCountAna(resultats,choix="dist")
		plot.WordCountAna(resultats,choix="cons")
	}
	return(resultats)
}
