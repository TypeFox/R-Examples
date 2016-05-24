######
check.distmat=function(obj,tol=1e-9) {
	if(class(obj)=="matrix" & length(unique(dim(obj)))==1) {
		if(all(diag(obj)-0<=tol)) {
			return(TRUE)
		} else {
			return(FALSE)
		}
	} else {
		return(FALSE)
	}
}


######
match.spacodi.data<-function(sp.plot, phy=NULL, sp.traits=NULL, prune=TRUE, verbose=FALSE) {
# major error checking
	sp.plot=as.matrix(sp.plot)
	if(is.null(row.names(sp.plot))) stop("Check that sp.plot has row names as species.") 
	if(is.null(colnames(sp.plot))) {
		warning("sp.plot does not appear to have plots as column names.")
		names(sp.plot)=paste("plot",1:ncol(sp.plot),sep="")
	}
	
	if(!missing(phy)) {
		if(class(phy)=="phylo") { 
			if(length(unique(phy$tip.label))!=length(phy$tip.label)) stop("Redundant taxa were found in tree.")
			if(!is.null(phy$node.label)) phy$node.label=NULL
		}
	}
	
	if(length(unique(names(sp.plot)))!=length(names(sp.plot))) stop("Redundant plots were found in sp.plot.")
	
# check poor values in sp.plot
	if(any(!is.finite(sp.plot))) {
		poor=which(!is.finite(sp.plot))
		poor%%nrow(sp.plot)->poor.rows
		ceiling(poor/nrow(sp.plot))->poor.cols
		poor.data=as.list(paste(rownames(sp.plot)[poor.rows],colnames(sp.plot)[poor.cols]))
		stop("Poor data values found in supplied sp.plot (NA, NaN, or Inf):\n\n\t",toString(poor.data),"\n")
	}
	
# find undersampled plots
	prune.sp=function(sp.plot, verbose=FALSE){
		rr=rownames(sp.plot)
		l.spp=nrow(sp.plot)
		drop.plots=vector()
		for(sp in 1:ncol(sp.plot)) {
			l.nulls=length(which(sp.plot[,sp]==0))
			if((l.spp-l.nulls)<2) {
				drop.plots=cbind(drop.plots, sp)
			}
		}
		
		plot.names.orig=colnames(sp.plot)
		dropped.plots=plot.names.orig[drop.plots]
		if(length(drop.plots)!=0) {
			if(verbose)message({cat("\nThe following plots were dropped from sp.plot:\n\t");cat(dropped.plots, sep=" "); cat("\n")})
			sp.plot=as.matrix(sp.plot[,-as.numeric(drop.plots[!is.na(drop.plots)])])
			rownames(sp.plot)=rr
			colnames(sp.plot)=plot.names.orig[which(!plot.names.orig%in%plot.names.orig[drop.plots])]
			if(is.null(ncol(sp.plot))) sp.plot=NULL
		}
		return(sp.plot)
	}
	
	if(prune) sp.plot=prune.sp(sp.plot, verbose=verbose)
	
# match ordering and size of all data objects	
	usable.species=rownames(sp.plot)
	
# phy only
	if(missing(sp.traits) & !missing(phy)) {
		if(class(phy)=="phylo") {
			usable.species.tmp=intersect(phy$tip.label, usable.species)
			phy=reorderspacodiobj(phy,usable.species.tmp)
			usable.species=usable.species.tmp[match(phy$tip.label, usable.species.tmp)]
			if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
			r.out=list(sp.plot=sp.plot,sp.tree=phy)
		} else if(check.distmat(phy)) {
			usable.species.tmp=intersect(rownames(phy), usable.species)
			phy=reorderspacodiobj(phy,usable.species.tmp)
			usable.species=usable.species.tmp[match(rownames(phy), usable.species.tmp)]
			if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
			r.out=list(sp.plot=sp.plot,sp.tree=phy)			
		}
	}
	
# sp.traits only
	if(!missing(sp.traits) & missing(phy)) {
		usable.species=intersect(usable.species, rownames(sp.traits))
		if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
		r.out=list(sp.plot=sp.plot,sp.traits=reorderspacodiobj(sp.traits, usable.species))			
	}
	
# both sp.traits and phy
	if(!missing(sp.traits) & !missing(phy)) {
		if(class(phy)=="phylo") {			
			usable.species.tmp=intersect(phy$tip.label, intersect(rownames(sp.traits), usable.species))
			phy=reorderspacodiobj(phy,usable.species.tmp)
			usable.species=usable.species.tmp[match(phy$tip.label, usable.species.tmp)]
			if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
			r.out=list(sp.plot=sp.plot,sp.tree=reorderspacodiobj(phy, usable.species), sp.traits=reorderspacodiobj(sp.traits, usable.species))
		} else if(check.distmat(phy)) {
			usable.species.tmp=intersect(rownames(phy), intersect(rownames(sp.traits), usable.species))
			phy=reorderspacodiobj(phy,usable.species.tmp)
			usable.species=usable.species.tmp[match(rownames(phy), usable.species.tmp)]
			if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
			r.out=list(sp.plot=sp.plot,sp.tree=reorderspacodiobj(phy, usable.species), sp.traits=reorderspacodiobj(sp.traits, usable.species))			
		}
	}
	
# neither sp.traits nor phy
	if(missing(sp.traits) & missing(phy)) {
		if(prune) sp.plot=prune.sp(reorderspacodiobj(sp.plot, usable.species),verbose=verbose) else sp.plot=reorderspacodiobj(sp.plot, usable.species)
		r.out=list(sp.plot=sp.plot)
	}
	
	return(r.out)
}	


######
reorderspacodiobj=function(obj, names) {
	if(class(obj)=="phylo") {
		obj.labels=obj$tip.label
		if(any(!obj.labels%in%names)) obj=drop.tip(obj, obj.labels[!obj.labels%in%names]) else obj=obj
	} else if(check.distmat(obj)) {
		obj.labels=rownames(obj)
		obj=as.matrix(obj[match(names,obj.labels),match(names,obj.labels)])
		rownames(obj)<-colnames(obj)<-names
	} else {
		nn=colnames(obj)
		obj.labels=rownames(obj)
		obj=as.data.frame(obj[match(names,obj.labels),])
		names(obj)=nn
	}
	if(!all(names%in%obj.labels)) warning(paste(paste(names[!names%in%obj.labels],sep=" ",collapse=" "),"were not found in the supplied object",sep=" ")) 
	return(obj)
}


######
write.spacodi.data <-
function(sp.plot, outfile){
	if(file.exists(outfile)){
		warning("Overwrote existing outfile.")
		unlink(outfile)
	}
	names=names(sp.plot)
	for(n in 1:length(names)){cat(c("\t",names[n]),file=outfile,append=TRUE,sep="")}
	cat("\n",file=outfile,append=TRUE,sep="")
	write.table(sp.plot,outfile,quote=FALSE,col.names=FALSE,append=TRUE,sep="\t")
}


######
as.phylocom <-
function(data, picante=FALSE, outfile=NULL){
	if(picante) data=as.spacodi(data)
	
	if(class(data)!="data.frame"){
		if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
	}
	if(length(unique(row.names(data)))!=nrow(data))warning("Data do not appear to be in proper format.  Results may be nonsensical.")
	plots=ncol(data)
	spp=nrow(data)
	out=array(dim=c(spp*plots,3))
	for(plot in 1:plots){
		start=((plot-1)*spp)+1
		end=plot*spp
		out[start:end,]=cbind(names(data)[plot], data[,plot], row.names(data))
	}
	out=as.data.frame(out)
	names(out)=c("plot","samples","species")
	if(!is.null(outfile)){write.table(out,outfile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")}
	return(out)
}


######
as.picante <-
function(data, outfile=NULL){
	if(ncol(data)==3) {
		message("Formatting assumed to be that used with phylocom")
		if(class(data)!="data.frame"){
			if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
		}
		species=as.character(unique(data[,3]))
		dd=split(data,data[,1])
		out=array(dim=c(length(species)->spp,length(dd)->plots))
		for(plot in 1:plots){
			cur.array=array(dim=c(spp, 1))
			ind=as.numeric(as.vector(dd[[plot]][,2]))
			dd[[plot]]->cur.plot
			names(ind)=as.character(cur.plot[,3])
			for(r in 1:nrow(cur.plot)){
				cur.array[which(names(ind[r])==species)]=ind[r]
			}
			cur.array[which(is.na(cur.array))]=0
			out[,plot]=as.numeric(cur.array)
		}
		out=as.data.frame(t(out))
		names(out)=species
		row.names(out)=names(dd)
		if(!is.null(outfile)){write.spacodi.data(out,outfile)}
		return(out)
	} else {
		message("Formatting assumed to be that used with spacodi")
		if(class(data)!="data.frame"){
			if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
		}
		if(is.null(row.names(data)))row.names(data)=paste("plot",seq(1:nrow(data)))
		if(is.null(names(data)))names(data)=paste("species",seq(1:nrow(data)))
		
		out=as.matrix(t(data))
		if(!is.null(outfile)){write.table(out,outfile,quote=FALSE)}
		return(out)
	}
}


######
as.spacodi <-
function(data, outfile=NULL){
	if(ncol(data)!=3) {
		message("Formatting assumed to be that used with picante")
		if(class(data)!="data.frame"){
			if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
		}
		
		if(is.null(row.names(data)))row.names(data)=paste("plot",seq(1:nrow(data)))
		if(is.null(names(data)))names(data)=paste("species",seq(1:nrow(data)))
		
		out=as.matrix(t(data))
		if(!is.null(outfile)){write.table(out,outfile,quote=FALSE)}
		return(out)		
	} else {
		message("Formatting assumed to be that used with phylocom")
		if(class(data)!="data.frame"){
			if(class(data)=="matrix") data=as.data.frame(data) else if(class(data)=="vector") data=as.data.frame(data) else stop("Function cannot handle data that are not in data.frame format")
		}
		species=as.character(unique(data[,3]))
		dd=split(data,data[,1])
		out=array(dim=c(length(species)->spp,length(dd)->plots))
		for(plot in 1:plots){
			cur.array=array(dim=c(spp, 1))
			ind=as.numeric(as.vector(dd[[plot]][,2]))
			dd[[plot]]->cur.plot
			names(ind)=as.character(cur.plot[,3])
			for(r in 1:nrow(cur.plot)){
				cur.array[which(names(ind[r])==species)]=ind[r]
			}
			cur.array[which(is.na(cur.array))]=0
			out[,plot]=as.numeric(cur.array)
		}
		out=as.data.frame(out)
		names(out)=names(dd)
		row.names(out)=species
		if(!is.null(outfile)){write.spacodi.data(out,outfile)}
		return(out)
	}
}

