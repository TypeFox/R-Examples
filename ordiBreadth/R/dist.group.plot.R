dist.group.plot <-
function(specialization,id,cex=1,PCoA=c(1,2),seg.col="pink",seg.wd=2,seg.lty=1,
pt.col="red",pt.pch=19,pt.cex=1.5,x.lim=NULL,y.lim=NULL,plant.plot="all.names",
rel.pch=19,rel.cex=1.5,rel.col="red",nrel.pch=19,nrel.cex=1.5,nrel.col="red",
verbose=TRUE,scaled=TRUE){
	
	plot(specialization$plants.ord[,PCoA[1]],specialization$plants.ord[,PCoA[2]],type="n",las=1,xlab=paste("PCoA ",PCoA[1]),ylab=paste("PCoA ",PCoA[2]),main="",ylim=y.lim,xlim=x.lim)
	if(verbose==TRUE){		
	mtext(paste("species",id,specialization$species[id]),line=2)
	if(scaled==TRUE){mtext(paste("#host plants:",sum(specialization$group.vectors[id,]==TRUE),"   scaled host breadth:",round(specialization$sca[id],digits=3)))}else{mtext(paste("#host plants:",sum(specialization$group.vectors[id,]==TRUE),"   host breadth:",round(specialization$tot[id],digits=3)))}
	}
	segs<-which(specialization$group.vectors[id,]==TRUE)
		for(i in segs){
			if(is.na(specialization$centroids.group[id,1])==TRUE){break}
			segments(specialization$plants.ord[i,PCoA[1]],specialization$plants.ord[i,PCoA[2]],specialization$centroids.group[id,PCoA[1]],specialization$centroids.group[id,PCoA[2]],col=seg.col,lwd=seg.wd,lty=seg.lty)
			}
		if(is.na(specialization$centroids.group[id,1])==TRUE){
				points(specialization$plants.ord[which(specialization$group.vectors[id,]==TRUE),PCoA[1]],specialization$plants.ord[which(specialization$group.vectors[id,]==TRUE),PCoA[2]],pch=pt.pch,cex=pt.cex,col=pt.col)}else{
				points(specialization$centroids.group[id,PCoA[1]],specialization$centroids.group[id,PCoA[2]],pch=pt.pch,cex=pt.cex,col=pt.col)
			}
	if(plant.plot=="all.names")text(specialization$plants.ord[,PCoA[1]],specialization$plants.ord[,PCoA[2]],rownames(specialization$plants.ord),cex=cex)
	if(plant.plot=="relevant")text(specialization$plants.ord[segs,PCoA[1]],specialization$plants.ord[segs,PCoA[2]],rownames(specialization$plants.ord[segs,]),cex=cex)
	if(plant.plot=="points"){
		points(specialization$plants.ord[segs,PCoA[1]],specialization$plants.ord[segs,PCoA[2]],pch=rel.pch,cex=rel.cex,col=rel.col)
		points(specialization$plants.ord[-segs,PCoA[1]],specialization$plants.ord[-segs,PCoA[2]],pch=nrel.pch,cex=nrel.cex,col=nrel.col)
	}
	}
