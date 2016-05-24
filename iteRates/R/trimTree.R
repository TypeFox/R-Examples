trimTree <-
function(phy,Time){
  	if (inherits(phy, "phylo")==FALSE){stop("object is not class phylo\n")}
    	if(is.ultrametric(phy)==FALSE){stop("\nTree not ultrametric\n")}
	if(max(branching.times(phy))<=Time){stop("\nTruncation time exceeds age of tree\n")}
	phy<-multi2di(phy,random=FALSE)
	phy<-read.tree(text=write.tree(phy))
		trees<-list()
		trees[[1]]<-phy
		newlabs<-trees[[1]]$tip.label
		j <- 2
	repeat{
		phy$tip.label<-1:length(phy$tip.label)
		phy$node.label<-(length(phy$tip.label)+1):(length(phy$tip.label)+length(phy$tip.label)-1)###
		rowstips<-which(phy$edge[,2]<=length(phy$tip.label))
		Realedge<-phy$edge.length
	
			if(sum((phy$edge.length[rowstips]-Time)<0)>0){
				phy$edge.length[rowstips]<-phy$edge.length[rowstips]-Time
				rowClip<-which(phy$edge.length<0)
				rowclipnod<-unique(phy$edge[rowClip,1])
				ns<-array(dim=c(length(rowclipnod),2))
					for (i in 1:length(rowclipnod)){
						ns[i,]<-node.sons(phy,rowclipnod[i])
						}
				todrop<-which(apply(ns,1,diff)==1)
				phy$edge.length<-Realedge
				phy<-drop.tip2.6(phy,tip=ns[todrop,1])
				phy$tip.label<-1:length(phy$tip.label)
				phy$node.label<-(length(phy$tip.label)+1):(length(phy$tip.label)+length(phy$tip.label)-1)
				
					#collapse taxa names
					x<-paste(newlabs[ns[todrop,1]],newlabs[ns[todrop,2]],sep=" & ")
					newlabs[ns[todrop,1]]<-NA
					newlabs[ns[todrop,2]]<-x
					newlabs<-as.vector(na.omit(newlabs))
				
				trees[[j]]<-phy
				j=j+1
			}
			else{			
				phy$edge.length<-Realedge
				break}
		}
	tphy<-trees[[length(trees)]]
	tphy$tip.label<-1:length(tphy$tip.label)
	tphy$node.label<-(length(tphy$tip.label)+1):(length(tphy$tip.label)+length(tphy$tip.label)-1)###
	rowstips<-which(tphy$edge[,2]<=length(tphy$tip.label))
	tphy$edge.length[rowstips]<-tphy$edge.length[rowstips]-Time
	tphy$Ntip<-length(tphy$tip.label)
	tphy$new.tip.clades<-newlabs
	return(list(o.tree=trees[[1]],t.tree=tphy))	
	}

