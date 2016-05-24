######
resamp.1a <-
function(obj, abund.class.ratio = 4) {
	if(abund.class.ratio<=1)stop("Supplied abundance-class ratio did not appear sensible: choose a value greater than 1.")
	orig=obj
	abund <- rowSums(obj)
	n.spp   <- length(abund)
	
	aa <- abund.class.ratio
	
	while(1) {
		classes <- aa^(0:ifelse(aa<2, n.spp*(1/(aa-1)), n.spp))
		if(length(classes)>=n.spp) {
			classes <- unique(round(runif(1)*classes))
			break()
		}
	}
	
	class <- rep(NA, n.spp)
	for(i in 1:length(classes)){
		class[abund > classes[i]] <- i
	}
	if(any(is.na(class))) class[which(is.na(class))]=1
	new_name <- rep(NA, n.spp) 
	for(i in unique(class)){
		new_name[class == i] <- sample(rownames(obj[class == i,]))
	}
	row.names(obj) <- new_name
	return(obj[order(match(row.names(obj),row.names(orig))),])
}


######
resamp.1s <-
function(obj) {
	orig=obj
	row.names(obj) <- sample(row.names(obj))
	return(obj[match(row.names(orig),row.names(obj)),])
}


######
resamp.2s <-
function(obj) {
	for(nn in 1:ncol(obj)){
		obj[,nn]=sample(obj[,nn])
	}
	return(obj)
}


######
.gswapcheck=function(mat){
	if(sum(cc<-mat==0)==2){
		d=diag(cc)
		if(all(d) | !any(d)) return(TRUE) else return(FALSE)
	} else if(any(cc)){
		if(all(cc)) return(TRUE) else return(FALSE)
	} else {
		return(TRUE)
	}
}

######
resamp.2x=function (obj, level = 0.1)
{
    swaps = max(round(level * ncol(obj) * nrow(obj)),1)
    X = sum(obj)
    orig=obj
    for (swap in 1:swaps) {
    	while(1){
    		rcol = sample(1:ncol(obj), 2)
        	rspp = sample(1:nrow(obj), 2)
        	if(.gswapcheck(obj[rspp,rcol])){
        		obj[rspp[1], rcol[1]] = orig[rspp[2], rcol[1]]
        		obj[rspp[2], rcol[1]] = orig[rspp[1], rcol[1]]
        		obj[rspp[1], rcol[2]] = orig[rspp[2], rcol[2]]
        		obj[rspp[2], rcol[2]] = orig[rspp[1], rcol[2]]
        		orig = obj
    			break()
        	}
     	}
    }
    if (sum(obj) != X)
    warning("A poor result is likely.")
    return(obj)
}


######
resamp.3i <-
function(obj) {
	for(ss in 1:nrow(obj)){
		obj[ss,]=sample(obj[ss,])
	}
	return(obj)
}


######
resamp.3t <-
function(obj, dmat=NULL) {
	if(is.null(dmat))flag=TRUE else flag=FALSE
	names.orig=names(obj)
	names(obj)=seq(1:ncol(obj))
	if(is.null(dmat)) {
		dmat=as.data.frame(matrix(0,ncol(obj),ncol(obj)))
		names(dmat)=names.orig
		row.names(dmat)=names.orig
	}
	dmat=as.data.frame(as.matrix(dmat))
	if(!all(names(dmat)%in%names.orig) || ncol(obj)!=ncol(dmat) || ncol(dmat)!=nrow(dmat))stop("Names in distance matrix do not correspond to plot names")
	row.names(dmat)=names(obj)
	names(dmat)=names(obj)
	
# find all distances from plot.tt to plot.tt + some shifter value (e.g., '3' would be plot1 to plot4, plot10 to plot3, ... plotN to plotN+3)
# tabulate these values and find the average distance from each plot to every plot+'shifter'
	torus=rep(1:ncol(obj),2)
	plus.array=array(dim=c(1,ncol(obj)))
	torus.array=array(dim=c(ncol(obj), ncol(obj)))
	for(plus in 1:ncol(obj)) {
		for(tt in 1:ncol(torus.array)) {
			from=tt
			to=torus[tt+plus]
			d.tt=dmat[from, to]
			torus.array[tt,plus]=d.tt
		}
		plus.array[1,plus]=mean(torus.array[,plus],na.rm=TRUE)
	}
	
	plus.array=as.data.frame(plus.array)
	names(plus.array)=names(obj)
	plus.array=plus.array[,order(plus.array)]
	
# randomly generate a value between 0 and maximum average distance from a plot to every other plot+'shifter'
# shift species abundances by the randomly chosen 'shifter' 
	for(ss in 1:nrow(obj)){
		if(!flag) {
			shifter=as.numeric(names(plus.array))[min(which(plus.array>=runif(1,min=min(plus.array), max=max(plus.array))))]
		} else {shifter = sample(as.numeric(names(plus.array)),1)}
		
		t.array=array(dim=c(ncol(obj),2))
		t.array[,1]=1:ncol(obj)
		for(o in 1:ncol(obj)){
			tt=torus[shifter+(o-1)]
			t.array[tt,2]=obj[ss,o]
		}
		obj[ss,]=t.array[,2]
	}
	if(flag) message("Plots were assumed to be equidistant from one another.")
	res=obj
	names(res)=names.orig
	return(res)
}

######
resamp.3x=function (obj, level = 0.1)
{
    swaps = max(round(level * ncol(obj) * nrow(obj)),1)
    X = sum(obj)
    orig=obj
    for (swap in 1:swaps) {
    	while(1){
    		rcol = sample(1:ncol(obj), 2)
        	rspp = sample(1:nrow(obj), 2)
        	if(.gswapcheck(obj[rspp,rcol])){
        		obj[rspp[1], rcol[1]] = orig[rspp[1], rcol[2]]
        		obj[rspp[1], rcol[2]] = orig[rspp[1], rcol[1]]
        		obj[rspp[2], rcol[1]] = orig[rspp[2], rcol[2]]
        		obj[rspp[2], rcol[2]] = orig[rspp[2], rcol[1]]
        		orig = obj
        		break()
        	}
     	}
    }
    if (sum(obj) != X)
    warning("A poor result is likely.")
    return(obj)
}


######
resamp.phy <-
function(phy, node=NULL, time.threshold=1, proportion=TRUE) {
	
	if(all(is.numeric(node))) {
	# allow tip reshuffling by prespecified node
		
		if(all(node>Ntip(phy)) & all(node<=max(phy$edge))) {
			internals=node
		} else {
			stop("Supplied 'node' value(s) appear misspecified.") 
		}
	} else {
	# reshuffling of tips by temporal constraint

		n=Ntip(phy)
		b=branching.times(phy)
		if(proportion) b=b/max(b)
		tt=sapply(b,function(x) withinrange(x,0,time.threshold))
		if(any(tt)){
			tt=as.numeric(names(tt)[which(tt)])
			affected.nodes=c()
			internals=c()
			to.do=lapply(tt, function(x) get.descendants.of.node(x, phy, tips=FALSE))
			for(nn in 1:length(to.do)) {
				if(any(!to.do[[nn]]%in%affected.nodes)) {
					internals=c(internals, tt[nn])
				}
				affected.nodes=c(affected.nodes, to.do[[nn]][!to.do[[nn]]%in%affected.nodes])
			}
		} else {
			warning("No resampling possible")
			return(NULL)
		}	
	}
	
	# SHUFFLE TIPS for NON-NESTED CLADES
	to.rand=lapply(internals, function(x) get.descendants.of.node(x, phy, tips=TRUE))
	for(nn in 1:length(to.rand)) {
		phy$tip.label[to.rand[[nn]]]=sample(phy$tip.label[to.rand[[nn]]])
	}
	
	return(phy)
}


