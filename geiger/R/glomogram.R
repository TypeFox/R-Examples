
.replace.tip.with.subtree=function(phy, subtree, tip, split_age, tol=.Machine$double.eps^0.5){
	
	stem=which(phy$tip.label==tip)
	
	if(!length(stem)) {
		print(tip)
		print(subtree)
		stop("Cannot find 'tip' in 'phy'")
	}
	if(!stem<=Ntip(phy)) stop("Replaced 'stem' must be a tip.")
	if(!class(subtree)=="phylo") stop("'subtree' must be a 'phylo' object.")
	if(!tip%in%phy$tip.label) stop("'tip' must occur in 'phy$tip.label'.")
#	if(!is.ultrametric(phy, tol=tol)) stop("'phy' must be ultrametric.")
	if(any(subtree$tip.label%in%phy$tip.label)) stop("Found matching tips in 'subtree' and 'phy'.")
	
	phy=reorder(phy)
	subtree=reorder(subtree)
	
	get.edgelength=function(node, phy) if(node==Ntip(phy)+1) return(0) else return(phy$edge.length[phy$edge[,2]==node])
	
	stem_length=get.edgelength(stem, phy)
	a=.get.ancestor.of.node(stem, phy)
	a_length=get.edgelength(a, phy)
	
	#	if(split_age>(stem_length+a_length)) stop("'split_age' is inconsistent with edge lengths in 'phy'")
	## FIXME: allow subtree to be pasted more rootward than node time of 'stem'? 
	
	if(split_age>(stem_length)) stop("'split_age' is inconsistent with edge lengths in 'phy'")
	
	# rework edge matrices
	n=Ntip(phy)
	sn=Ntip(subtree)
	N=n+sn
	a=a+sn
	edge=phy$edge
	es.t<- edge>stem & edge<=n
	edge[es.t]=edge[es.t]-1
	lengths=phy$edge.length
	sedge=subtree$edge
	en.t<- edge>n
	edge[en.t]=edge[en.t]+sn
	ea.t<- edge>a
	edge[ea.t]=edge[ea.t]+sn-1
	ss.t<- sedge>sn
	sedge[ss.t]=sedge[ss.t]+a-sn
	sn.t<- sedge<=sn
	sedge[sn.t]=sedge[sn.t]+n-1
	
	insert_edge=sedge
	insert_length=subtree$edge.length
	
	stem_slot=which(phy$edge[,2]==stem)
	
	graft_replacement.phy=function(edge, edge.length, phy, insert_edge, insert_length, subtree, slot, stem){
		
		newedge.tmp=rbind(edge[1:slot,], insert_edge)
		newlengths.tmp=c(edge.length[1:slot], insert_length)
		if(slot==nrow(edge)){
			newedge=newedge.tmp
			newlengths=newlengths.tmp
		} else {
			newedge=rbind(newedge.tmp, edge[(slot+1):nrow(edge),])
			newlengths=c(newlengths.tmp, edge.length[(slot+1):length(edge.length)])
		}
		
		newnnode=Nnode(phy)+Nnode(subtree)
		
		newphy=list(edge=newedge, edge.length=newlengths, tip.label=c(phy$tip.label[-stem], subtree$tip.label), Nnode=newnnode)
		newphy
	}
	
	newphy=graft_replacement.phy(edge, lengths, phy, insert_edge, insert_length, subtree, stem_slot, stem)
	newphy$edge[stem_slot,2]=a+1
	newphy$edge.length[stem_slot]=newphy$edge.length[stem_slot]-split_age
	nl.t<- newphy$edge>length(newphy$tip.label)
	newphy$edge[nl.t]=newphy$edge[nl.t]-1
	class(newphy)="phylo"
	newphy=reorder(newphy)
	return(newphy)
}

.prep.glomogram=function(labels, phy){
	N=Ntip(phy)
	dd=.cache.descendants(phy)
	nn=c(phy$tip.label, phy$node.label)
	mm=match(labels, nn)
	if(any(is.na(mm))) stop("Some 'labels' not encountered in 'phy'")
	tmp=dd$adesc[mm]
	ck=sapply(1:length(tmp), function(idx) {
		   cur=tmp[[idx]]
		   any(cur%in%mm[-idx])
		   })
	if(any(ck)) stop("Nested 'labels' encountered")
	names(mm)=labels
	inner=mm[mm>N]
	if(length(inner)){
		tips=dd$tips[inner]
		for(i in 1:length(inner)){
			phy$tip.label[tips[[i]]]=names(inner)[i]
		}		
	}
	.unique.phylo(phy)
}


glomogram.phylo=function(phy, subtrees){
# subtrees: a named list of subtrees (each can be a 'multiPhylo' object)
#		-- minimally has 'subtree' and 'age' elements
#		-- if 'tip' is missing, the name of the 'subtrees' object is taken as 'tip'
# phy: the skeleton tree with tip labels corresponding to names in 'subtrees'
	
	if(is.null(names(subtrees))) stop("'subtrees' must have associated names")
	
	phy=.prep.glomogram(names(subtrees), phy)
	
	if(class(phy)=="multiPhylo") phy=phy[[sample(1:length(phy),1)]]
	phy$node.label=NULL
	done=c()
	for(i in 1:length(subtrees)){
#		cat(paste(names(subtrees)[i], "\n", sep=""))
		
		cur=subtrees[[i]]
		
		if(class(cur)=="multiPhylo") {
			cur=cur[[sample(1:length(cur), 1)]]
		} 
		
		if(Ntip(cur)==1){
			phy$tip.label[phy$tip.label==cur$tip]=cur$tip.label
		} else {
			ht=heights(cur)
			curage=max(ht)
			curtip=names(subtrees)[i]
			if(!curtip%in%done){
				done=c(done, curtip)
				phy=.replace.tip.with.subtree(phy, cur, curtip, curage)		
			} else {
				warning("Duplicate tips encountered: current 'subtree' will be skipped.")
				next()
			}
		}
	}	
	
## ENSURE proper node IDs (problematic when dealing with polytomous subtrees)
## VERY SLOW -- FIXME
	xx=unique(sort(phy$edge))
	nn=1:length(xx)
	N=Ntip(phy)
	inner=xx[xx>N]
	for(i in 1:length(inner)){
		nd=inner[i]
		phy$edge[phy$edge==nd]=i+N
	}
	phy
}
