#' Convert Simulated Taxon Data into a Phylogeny
#' 
#' Converts temporal and ancestor-descendant relationships of taxa into a
#' time-scaled phylogeny
#' 
#' @details 
#' As described in the documentation for taxa2cladogram, the relationships
#' among morphotaxa in the fossil record are difficult to describe in terms of
#' traditional phylogenies. One possibility is to arbitrarily choose particular
#' instantaneous points of time in the range of some taxa and describe the
#' time-scaled relationships of the populations present at those dates. This is
#' the tactic used by taxa2phylo.
#' 
#' By default, the dates selected ('obs-time' argument) are the last occurances
#' of the taxon, so a simple use of this function will produce a time-scaled
#' tree which describes the relaitonships of the populations present at the
#' last occurance of each taxon in the sampled data. Alternatively, obs_time
#' can be supplied with different dates within the taxon ranges.
#' 
#' All data relating to when static morpho-taxa appear or disappear in the
#' record is lost; branching points will be the actual time of speciation,
#' which (under budding) will often be in the middle of the temporal range of a
#' taxon.
#' 
#' Cryptic taxa are not dropped or merged as can be done with taxa2cladogram.
#' The purpose of taxa2phylo is to obtain the 'true' pattern of evolution for
#' the observation times, independent of what we might actually be able to
#' recover, for the purpose of comparing in simulation analyses.
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' @param taxad A five-column matrix of taxonomic data, as output by
#' \code{fossilRecord2fossilTaxa} via simulations produced using \code{simFossilRecord}

#' @param obs_time A vector of per-taxon times of observation which must be in
#' the same order of taxa as in the object taxad; if NULL, the LADs (column 4)
#' in taxad2 are used

#' @param plot Plot the resulting phylogeny?

#' @return The resulting phylogeny with branch lengths is output as an object
#' of class phylo. This function will output trees with the element $root.time,
#' which is the time of the root divergence in absolute time.
#' 
#' The tip labels are the rownames from the simulation input; see documentation
#' for \code{simFossilRecord} and \code{fossilRecord2fossilTaxa} documentation for details.

#' @note DO NOT use this function to time-scale a real tree for a real dataset.
#' It assumes you know the divergence/speciation times of the branching nodes
#' and relationships perfectly, which is almost impossible given the
#' undersampled nature of the fossil record. Use timePaleoPhy or
#' cal3TimePaleoPhy instead.
#' 
#' DO use this function when doing simulations and you want to make a tree of
#' the 'true' history, such as for simulating trait evolution along
#' phylogenetic branches.
#' 
#' Unlike \code{taxa2cladogram}, this function does not merge cryptic taxa in output
#' from \code{simFossilRecord} (via \code{fossilRecord2fossilTaxa})
#' and I do not offer an option to secondarily drop them.
#' The tip labels should provide the necessary information for users to drop
#' such taxa, however. See \link{simFossilRecord}.

#' @author David W. Bapst

#' @seealso \code{\link{simFossilRecord}}, \code{\link{taxa2cladogram}}, \link{fossilRecord2fossilTaxa}

#' @examples
#' 
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'  	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
#' tree<-taxa2phylo(taxa)
#' phyloDiv(tree)
#' 
#' #now a phylogeny with tips placed at the apparent time of extinction for each taxon
#' rangesCont<-sampleRanges(taxa,r=0.5)
#' tree<-taxa2phylo(taxa,obs_time=rangesCont[,2])
#' phyloDiv(tree,drop.ZLB=FALSE)
#' #note that it drops taxa which were never sampled!
#' 
#' #testing with cryptic speciation
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, prop.cryptic=0.5, nruns=1,
#'  	nTotalTaxa=c(30,40), nExtant=0, count.cryptic=TRUE)
#' taxaCrypt<-fossilRecord2fossilTaxa(record)
#' treeCrypt<-taxa2phylo(taxaCrypt)
#' layout(1)
#' plot(treeCrypt)
#' axisPhylo()
#' 
#' @export taxa2phylo
taxa2phylo<-function(taxad,obs_time=NULL,plot=FALSE){
		#obs_time=NULL; plot=TRUE
	#INPUT a taxad object and a vector of observation times for each species
		#if obs=NULL, LADs in taxad1 are used as observation times
		#all times must be in backwards format (zero is present)
	#root.time
		#ALL TREES ARE OUTPUT WITH ELEMENTs "$root.time"
		#this is the time of the root on the tree, which is important for comparing across trees
		#this must be calculated prior to adding anything to terminal branches
	#OUTPUT an ape phylo object with the tips at the times of observation
	#require(ape)
	#important checks 06-17-15
	#check that taxa are ordered by FAD
	if(!identical(order(-taxad[,3]),1:nrow(taxad))){
		#let's coerce the taxad so that it satisfies this
		taxadNew<-taxad[order(-taxad[,3]),]
		if(!is.null(obs_time)){
			#reorder obs_time too!
			obs_time<-obs_time[order(-taxad[,3])]
			}
		#reassign ancestor IDs
			#DON'T FORGET YOU NEED TO REORDER THEM
		newAnc<-sapply(taxad[order(-taxad[,3]),2],function(x) {
			if(is.na(x)){
				NA
			}else{
				which(sapply(taxadNew[,1],identical,x))
				}
			})
		if(is.vector(newAnc)){
			taxadNew[,2]<-newAnc
		}else{
			stop("ancestor IDs cannot be reassigned properly")
			}
		#reassign taxon IDs
		taxadNew[,1]<-1:nrow(taxad)
		taxad<-taxadNew
		#
		#check
		if(!identical(order(-taxad[,3]),1:nrow(taxad))){
			stop("Cannot coerce input table so taxa are ordered with FADs going from first to last")
			}
		#message(paste0("Input table must be ordered with FADs going from first to last",
		#	" \n coercing taxon order to first row"))
		#
		#stop("input table must be ordered with FADs going from first to last")
		}
	# check root ancestor is in first row
	if(!is.na(taxad[1,2])){
		#find the damn root
		isRoot<-which(is.na(taxad[,2]))
		if(length(isRoot)!=1){
			if(length(isRoot)>1){
				stop("Multiple taxa listed as an apparent root (ancestor is NA)")}
			if(length(isRoot)<1){
				stop("No taxa are listed as an apparent root (i.e. ancestor is NA)")}
		}else{
			stop("Root taxon (ancestor listed as NA) must be in row 1")
			#newRoot<-taxad[isRoot,,drop=FALSE]
			#newRoot[,1]<-1
			#taxadDrop<-taxad[-isRoot,]
			#taxadDrop[taxadDrop[,1]==1,]<-isRoot
			#taxad<-rbind(newRoot,taxadDrop)
			#taxad<-taxad[order(taxad[,1]),]
			#message(paste0("Root ancestor (ancestor listed as NA) is not in row 1,",
			#	" \n coercing root ancestor to first row"))
			}		
		}
	taxad1<-taxad[,1:4,drop=FALSE]
	#some checks
	if(!testParentChild(parentChild=taxad1[,2:1])){
		stop("input anc-desc relationships are inconsistent")}
	#
	#convert time so it runs from forward in time
	if(any((taxad1[,4]-taxad1[,3])<0)){
		taxad1[,3:4]<-max(taxad1[,3:4])-taxad1[,3:4]}
	if(any((taxad1[,4]-taxad1[,3])<0)){
		stop("Last occurrences appear to occur before first occurrences?")}
	if(any(table(taxad1[,1])>1)){
		stop(paste("Duplicated values of taxon ID's in input:",
			which(tabulate(taxad1[,1])>1),collapse=" "))}
	if(is.null(obs_time)){
		obs<-taxad1[,4]
	}else{
		obs<-max(taxad[,3:4])-obs_time
		#check if the times of observations are outside of original taxon ranges
		#nameMatch<-match(names(time.obs),rownames(taxR))
		#if(any(is.na(nameMatch))){stop("ERROR: names on time.obs and in candleRes don't match")}
		obsOutRange<-sapply(1:length(obs_time),function(x)
			if(is.na(obs_time[x])){
				FALSE
			}else{
				(obs_time[x]>taxad[x,3])|(obs_time[x]<taxad[x,4])
			})
		if(any(obsOutRange)){
			stop(paste0("ERROR: Given obs_time are outside of the original taxon ranges!",
				"\n If cryptic taxa, perhaps you forgot to set merge.cryptic=FALSE?"))}
		}
	if(nrow(taxad1)!=length(obs)){
		stop("Number of observations are not equal to number of lineages!")}
	#make observations as fake taxa, assuming that observations are WITHIN actual taxon ranges
	fake_taxa<-matrix(sapply((1:nrow(taxad1))[!is.na(obs)],function(x)
		 c(nrow(taxad1)+x,taxad1[x,1],obs[x],obs[x])),,4,byrow=TRUE)
	fake_taxa[,1]<-(1:nrow(fake_taxa))+nrow(taxad1)
	taxad2<-rbind(taxad1,fake_taxa)
	ntaxa<-nrow(taxad2)
	#MAKE IT INTO AN NODE/EDGE-BASED PHYLOGENY
	#find descendents of every taxon
	desc<-lapply(taxad2[,1],function(x)
		taxad2[c(FALSE,taxad2[-1,2]==x),1])	
	#get time of desc births
	births2<-lapply(desc,function(x)
		 if(length(x)>0){
			sapply(x,function(y) taxad2[y==taxad2[,1],3])
		})
	desc<-lapply(1:length(desc),function(x)	#MORE STUPIDLY COMPLICATED CODE
		if(length(desc[[x]])>0){
			desc[[x]][match(1:length(births2[[x]]),rank(births2[[x]],ties.method="random"))]
		}else{
			numeric()
			})
	#time of desc births
	births<-lapply(desc,function(x)
		 if(length(x)>0){sapply(x,function(y) taxad2[y==taxad2[,1],3])})	
	#make events list: first event is taxon birth, with that taxon as desc,
		# next is desc births, extinction not recorded but implied
	events<-lapply(1:ntaxa,function(x)
		if(length(desc[[x]])>0){
			c(taxad2[x,1],desc[[x]])
		}else{
			c(taxad2[x,1])
			})
	#times of events: taxon birth, desc births, extinction
	times<-lapply(1:ntaxa,function(x)
		if(length(births[[x]])>0){
			c(taxad2[x,3],births[[x]],taxad2[x,4])
		}else{
			c(taxad2[x,3],taxad2[x,4])
		})
	#labels for each lineage segment between event times
		# (these may as well represent the daughter node ID too)
		#use decimals to keep track of segments, set first segment as X.0
		#as long as a single taxon doesn't have millions of descendants,
			# in which case the IDs may stop being unique...
	nseg<-sapply(times,length)-1
	seg_labs<-lapply(1:ntaxa,function(x) taxad2[x,1]+(1:nseg[x]/(nseg[x]+1)))
	seg_labs<-lapply(seg_labs,function(x) c(floor(x[1]),x[-1]))
	#now, find the mother segment for each taxon
	taxa_anc<-c(0,sapply(2:ntaxa,function(x)
		unlist(seg_labs[taxad2[x,2]==taxad2[,1]])[
			which(unlist(events[taxad2[x,2]==taxad2[,1]])==taxad2[x,1])-1]))
	#now make list of all seg ids for all mom segs
	moms2<-lapply(seg_labs,function(x) x[-length(x)])	
	moms<-lapply(1:ntaxa,function(x) c(taxa_anc[x],unlist(moms2[[x]])))
	lengths<-lapply(times,diff)	#get lengths of segments 
	term<-sapply(unlist(seg_labs),function(x) !any(unlist(moms)==x))	#which are terminal?
	#make edge data.frame, with id, anc-id, length and logical indicating terminal branches
	edgeD<-data.frame(id=unlist(seg_labs),anc=unlist(moms),brlen=unlist(lengths),term=term)
	MRCA<-min(edgeD[sapply(edgeD$id,function(x) 1<sum(x==edgeD$anc)),1])
	#edgeD[edgeD[,2]==MRCA,]
	edgeD<-edgeD[-which(edgeD[,1]<=MRCA),]	#want to drop any extraneous root as well
	#add/leave only the terminals needed for the observations, remove others
	droppers<-which(as.logical(edgeD[,4]) & edgeD[,1]<(nrow(taxad1)+1))
	while(length(droppers)>0){
		edgeD<-edgeD[-droppers,]
		edgeD[,4]<-sapply(edgeD[,1],function(x) !any(edgeD[,2]==x))
		droppers<-which(as.logical(edgeD[,4]) & edgeD[,1]<(nrow(taxad1)+1))
		}
	edgeD[,4]<-sapply(edgeD[,1],function(x) !any(edgeD[,2]==x))
	#collapse single internal nodes
	ndesc<-sapply(edgeD[,1],function(x) sum(x==edgeD[,2]))	#ndesc from each node
	while(any(ndesc==1)){		#picks only internal branches with 1 desc
		#pick a single, if matrix, use first one
		epick<-edgeD[ndesc==1,]
		if(is.data.frame(epick)){epick<-epick[1,]}	
		edesc<-edgeD[edgeD$anc==epick$id,]
		#remove desc
		edgeD<-edgeD[-which(edgeD$id==edesc$id),]
		#replace picked edge w/new edge
		newe<-data.frame(id=edesc$id, anc=epick$anc,
			brlen=(epick$brlen+edesc$brlen), term=edesc$term)
		edgeD[edgeD$id==epick$id,]<-newe
		ndesc<-sapply(edgeD$id,function(x) sum(x==edgeD$anc))
		}
	ndesc<-sapply(edgeD$id,function(x) sum(x==edgeD$anc))
	#replace stupid decimal edge placeholders with clean numbers
	e_fix<-numeric()
	#the species in the original data
	e_fix[edgeD$term]<-sapply(edgeD[edgeD$term,1],function(x)
		 which(fake_taxa[,1]==x)) 	
	e_fix[!edgeD$term]<-sum(edgeD$term)+1+(1:sum(!edgeD$term))
	ea_fix<-sapply(edgeD$anc,function(x)
		 ifelse(x!=MRCA,e_fix[edgeD$id==x],sum(edgeD$term)+1))
	#NOW MAKE A TREE
	tlabs<-rownames(taxad1)[!is.na(obs)]
	edgf<-cbind(ea_fix,e_fix)
	colnames(edgf)<-NULL
	#check one more time before you make it a tree
	if(any(is.na(edgf))){
		stop("NAs introduced into edge matrix?")}
	if(!testParentChild(parentChild=edgf)){
		stop("produced edge matrix is inconsistent")}
	#Now really make it a tree
	tree1<-list(edge=edgf,tip.label=tlabs,edge.length=edgeD[,3],Nnode=length(unique(edgf[,1])))
	class(tree1)<-"phylo"						#ITS A TREE!
	#tree<-reorder(collapse.singles(tree1),"cladewise") 	#REORDER IT
	#if(!testEdgeMat(tree)){stop("Edge matrix has inconsistencies")}
	#tree<-read.tree(text=write.tree(tree))
	tree<-cleanNewPhylo(tree1)
	#now, root.time should be the time of the first obs PLUS the distance
		# from the earliest tip to the root
	first_obs_time<-max(taxad[,3:4])-min(obs,na.rm=TRUE)
	tree$root.time<-first_obs_time+min(node.depth.edgelength(tree)[1:Ntip(tree)])
	#make it so that root.time-max node dist from root must be below zero
	rootOffset<-tree$root.time-max(node.depth.edgelength(tree))
	if(rootOffset<0){
		tree$root.time<-tree$root.time-rootOffset
		}
	#plot
	if(plot){
		plot(ladderize(tree),show.tip.label=FALSE)
		axisPhylo()
		}
	#
	return(tree)
	}
