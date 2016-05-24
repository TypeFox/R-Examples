#' Creating a Taxon-Tree from Taxonomic Data Downloaded from the Paleobiology Database
#'
#' This function creates phylogeny-like object of type \code{phylo} from the taxonomic information
#' recorded in a taxonomy download from the PBDB for a given group. Two different algorithms are provided,
#' the default being based on parent-child taxon relationships, the other based on the nested Linnean hierarchy.

#' @details
#' This function should not be taken too seriously. Many groups in the Paleobiology Database have
#' out-of-date or very incomplete taxonomic information. This function is meant to help visualize
#' what information is present, and by use of time-scaling functions, allow us to visualize the intersection
#' of temporal and phylogenetic, mainly to look for incongruence due to either incorrect taxonomic placements,
#' erroneous occurrence data or both. 
#'
#' Note however that, contrary to common opinion among some paleontologists, taxon-trees may be just as useful for 
#' macroevolutionary studies as reconstructed phylogenies (Soul and Friedman, in press.).

#' @param data A table of taxonomic data collected from the Paleobiology Database, using the taxa list option
#' with show=phylo. Should work with versions 1.1-1.2 of the API, with either the 'pbdb' or 'com' vocab. However,
#' as 'accepted_name' is not available in API v1.1, the resulting tree will have a taxon's *original* name and not
#' any formally updated name.

#' @param rank The selected taxon rank; must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'.

#' @param cleanTree By default, the tree is run through a series of post-processing, including having singles collapsed,
#' nodes reordered and being written out as a Newick string and read back in, to ensure functionality with ape functions
#' and ape-derived functions. If FALSE, none of this post-processing is done and users should beware, as such trees can
#' lead to hard-crashes of R.

#' @param method Controls which algorithm is used for calculating the taxon-tree: either \code{method =} \code{"parentChild"}
#' (the default option) which converts the listed binary parent-child taxon relationships from the input PBDB data,
#' or \code{method = "Linnean"}, which converts a taxon-tree by creating a table of the Linnean
#' taxonomic assignments (family, order, etc), which are provided when
#' option 'show=phylo' is used in PBDB API calls.

#' @param tipSet This argument only impacts analyses where the argument
#' \code{method =} \code{"parentChild"} is also used. This \code{tipSet} controls
#' which taxa are selected as tip taxa for the
#' output tree. The default \code{tipSet =} \code{"nonParents"} selects all child taxa which
#' are not listed as parents in \code{parentChild}. Alternatively, \code{tipSet = "all"}
#' will add a tip to every internal node with the parent-taxon name encapsulated in
#' parentheses.

#' @param solveMissing Under \code{method =} \code{"parentChild"}, what should \code{makePBDBtaxonTree} do about
#' multiple 'floating' parent taxa, listed without their own parent taxon information in the input
#'  dataset under \code{data}? Each of these is essentially a separate root taxon, for a different set
#'  of parent-child relationships, and thus poses a problem as far as returning a single phylogeny is
#'  concerned. If \code{solveMissing = NULL} (the default), nothing is done and the operation halts with
#'  an error, reporting the identity of these taxa. Two alternative solutions are offered: first,
#'  \code{solveMissing =} \code{"mergeRoots"} will combine these disparate potential roots and link them to an
#'  artificially-constructed pseudo-root, which at least allows for visualization of the taxonomic
#'  structure in a limited dataset. Secondly, \code{solveMissing =} \code{"queryPBDB"} queries the Paleobiology
#'  Database repeatedly via the API for information on parent taxa of the 'floating' parents, and continues
#'  within a \code{while()} loop until only one such unassigned parent taxon remains. This latter option may
#'  talk a long time or never finish, depending on the linearity and taxonomic structures encountered in the
#'  PBDB taxonomic data; i.e. if someone a taxon was ultimately its own indirect child in some grand loop by
#'  mistake, then under this option \code{makePBDBtaxonTree} might never finish. In cases where taxonomy is
#'  bad due to weird and erroneous taxonomic assignments reported by the PBDB, this routine may search all
#'  the way back to a very ancient and deep taxon, such as the Eukaryota taxon.
#' Users should thus use \code{solveMissing =} \code{"queryPBDB"} only with caution.

#' @param APIversion Version of the Paleobiology Database API used by \code{makePBDBtaxonTree} when
#' \code{solveMissing =} \code{"queryPBDB"}. The current default is "1.1", which is the only option available
#' as of 05/05/2015. In the future, the improved API version "1.2" will be released on the public
#' PBDB server, which will become the new default for this function, but the option to return to "1.1"
#' behavior will be retained for .

# @param cleanDuplicate If TRUE (\emph{not} the default), duplicated taxa of a
# taxonomic rank \emph{not} selected by argument \code{rank}
# will be removed silently. Only duplicates of the taxonomic rank of interest
# will actually result in an error message.

#' @return
#' A phylogeny of class 'phylo', where each tip is a taxon of the given 'rank'. See additional details
#' regarding branch lengths can be found in the sub-algorithms used to create the taxon-tree by this function:
#' \code{\link{parentChild2taxonTree}} and \code{\link{taxonTable2taxonTree}}.
#'
#' Depending on the \code{method}
#' used, either the element \code{$parentChild} or \code{$taxonTable} is added to the list structure of
#' the output phylogeny object, which was used as input for one of the two algorithms mentioned above.
#'
#' Please note that when applied to output from the taxa option of the API version 1.1, the taxon names
#' returned are the \emph{original} taxon names as 'accepted_name' is not available in API v1.1, while
#' under API v1.2, the returned taxon names should be the most up-to-date formal names for those taxa.
#' Similar issues also effect the identification of parent taxa, as the accepted name of the
#' parent ID number is only provided in version 1.2 of the API.

#' @seealso
#' Two other functions in paleotree are used as sub-algorithms by \code{makePBDBtaxonTree}
#' to create the taxon-tree within this function,
#' and users should consult their manual pages for additional details:
#'
#' \code{\link{parentChild2taxonTree}} and \code{\link{taxonTable2taxonTree}}
#'
#' Other functions for manipulating PBDB data can be found at \code{\link{taxonSortPBDBocc}},
#' \code{\link{occData2timeList}}, and the example data at \code{\link{graptPBDB}}.

#' @author David W. Bapst

#' @references
#' Soul, L. C., and M. Friedman. In Press. Taxonomy and Phylogeny Can Yield
#' Comparable Results in Comparative Palaeontological Analyses. \emph{Systematic Biology} 
#' (\href{http://sysbio.oxfordjournals.org/content/early/2015/03/23/sysbio.syv015.abstract}{Link})
#' 

#' @examples
#' \dontrun{
#' 
#' easyGetPBDBtaxa<-function(taxon,show=c("phylo","img","app")){
#' 	#let's get some taxonomic data
#' 	taxaData<-read.csv(paste0("http://paleobiodb.org/",
#' 		"data1.1/taxa/list.txt?base_name=",taxon,
#' 		"&rel=all_children&show=",
#'		paste0(show,collapse=","),"&status=senior"),
#'		stringsAsFactors=FALSE)
#' 	return(taxaData)
#' 	}
#' 
#' #graptolites
#' graptData<-easyGetPBDBtaxa("Graptolithina")
#' graptTree<-makePBDBtaxonTree(graptData,"genus",
#' 	method="parentChild", solveMissing="queryPBDB")
#' #try Linnean
#' graptTree<-makePBDBtaxonTree(graptData,"genus",
#' 	method="Linnean")
#' plot(graptTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(graptTree$node.label,adj=c(0,1/2))
#' 
#' #conodonts
#' conoData<-easyGetPBDBtaxa("Conodonta")
#' conoTree<-makePBDBtaxonTree(conoData,"genus",
#' 	method="parentChild", solveMissing="queryPBDB")
#' plot(conoTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(conoTree$node.label,adj=c(0,1/2))
#' 
#' #asaphid trilobites
#' asaData<-easyGetPBDBtaxa("Asaphida")
#' asaTree<-makePBDBtaxonTree(asaData,"genus",
#' 	method="parentChild", solveMissing="queryPBDB")
#' plot(asaTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(asaTree$node.label,adj=c(0,1/2))
#' 
#' #Ornithischia
#' ornithData<-easyGetPBDBtaxa("Ornithischia")
#' ornithTree<-makePBDBtaxonTree(ornithData,"genus",
#' 	method="parentChild", solveMissing="queryPBDB")
#' #try Linnean
#' #need to drop repeated taxon first: Hylaeosaurus
#' ornithData<-ornithData[-(which(ornithData[,"taxon_name"]=="Hylaeosaurus")[1]),]
#' ornithTree<-makePBDBtaxonTree(ornithData,"genus",
#' 	method="Linnean")
#' plot(ornithTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(ornithTree$node.label,adj=c(0,1/2))
#' 
#' #Rhynchonellida
#' rynchData<-easyGetPBDBtaxa("Rhynchonellida")
#' rynchTree<-makePBDBtaxonTree(rynchData,"genus",
#' 	method="parentChild", solveMissing="queryPBDB")
#' plot(rynchTree,show.tip.label=FALSE,no.margin=TRUE,edge.width=0.35)
#' nodelabels(rynchTree$node.label,adj=c(0,1/2))
#' 
#' #some of these look pretty messy!
#' 
#' }
#' 
#' ###################################
#' \donttest{
#' 
#' #let's try time-scaling the graptolite tree
#' 
#' #get some example occurrence and taxonomic data
#' data(graptPBDB)
#' 
#' #get the taxon tree: Linnean method
#' graptTree<-makePBDBtaxonTree(graptTaxaPBDB, "genus", method="Linnean")
#' plot(graptTree,cex=0.4)
#' nodelabels(graptTree$node.label,cex=0.5)
#'
#' #get the taxon tree: parentChild method
#' graptTree<-makePBDBtaxonTree(graptTaxaPBDB, "genus", method="parentChild")
#' plot(graptTree,cex=0.4)
#' nodelabels(graptTree$node.label,cex=0.5)
#' 
#' #get time data from occurrences
#' graptOccGenus<-taxonSortPBDBocc(graptOccPBDB,rank="genus",onlyFormal=FALSE)
#' graptTimeGenus<-occData2timeList(occList=graptOccGenus)
#' 
#' #let's time-scale the parentChild tree with paleotree
#'		# use minimum branch length for visualization
#' 		# and nonstoch.bin so we plot maximal ranges
#' timeTree<-bin_timePaleoPhy(graptTree,timeList=graptTimeGenus,
#' 	nonstoch.bin=TRUE,type="mbl",vartime=3)
#' 
#' #drops a lot of taxa; some of this is due to mispellings, etc
#' 
#' }
#' \dontrun{
#' 
#' #make pretty plot with library strap
#' library(strap)
#' geoscalePhylo(timeTree, ages=timeTree$ranges.used)
#' nodelabels(timeTree$node.label,cex=0.5)
#'
#' }
#' 

#' @name makePBDBtaxonTree
#' @rdname makePBDBtaxonTree
#' @export
makePBDBtaxonTree<-function(data,rank,method="parentChild",solveMissing=NULL,
					tipSet="nonParents",cleanTree=TRUE,APIversion="1.1"){		
	# 
	# library(paleotree);data(graptPBDB);
	# data<-graptTaxaPBDB; rank="genus"; method="parentChild"; tipSet="nonParents"; cleanTree=TRUE; solveMissing=NULL
	# data<-graptTaxaPBDB; rank="genus"; method="parentChild"; tipSet="nonParents"; cleanTree=TRUE; solveMissing="queryPBDB"
	# data<-graptTaxaPBDB; rank="genus"; method="parentChild"; tipSet="nonParents"; cleanTree=TRUE; solveMissing="mergeRoots"
	# data<-graptTaxaPBDB; rank="genus"; method="Linnean"; 
	#
	#CHECKS
	if(length(method)!=1 | !is.character(method)){
		stop("method must be a single character value")}
	if(!any(method==c("Linnean","parentChild"))){
		stop("method must be one of either 'Linnean' or 'parentChild'")}
	if(!is.null(solveMissing)){
		if(length(solveMissing)>1 | !is.character(solveMissing)){
			stop("solveMissing must be either NULL or a single character value")}
		if(is.na(match(solveMissing,c("queryPBDB","mergeRoots")))){
			stop('solveMissing but be either NULL or "queryPBDB" or "mergeRoots"')}
		}
	if(!is(data,"data.frame")){stop("data isn't a data.frame")}
	if(length(rank)!=1 | !is.character(rank)){
		stop("rank must be a single character value")}
	if(!any(sapply(c("species","genus","family","order","class","phylum"),function(x) x==rank))){
		stop("rank must be one of 'species', 'genus', 'family', 'order', 'class' or 'phylum'")}
	#
	#translate to a common vocabulary
	data1<-translatePBDBtaxa(data)
	data1<-apply(data1,2,as.character)
	#
	if(method=="parentChild"){
		#need two things: a table of parent-child relationships as IDs
			#and a look-up table of IDs and taxon names
		#
		#filer out lower than selected rank
		# translate rank and taxon_rank to a number
		taxRankPBDB<-c("subspecies","species","subgenus","genus","subtribe","tribe","subfamily",
			"family","superfamily","infraorder","suborder","order","superorder","infraclass",
			"subclass","class","superclass","subphylum","phylum","superphylum","subkingdom",
			"kingdom","unranked clade","informal")	#keep informal as high, never drop!
		rank1<-which(rank==taxRankPBDB)
		numTaxonRank<-sapply(data1[,"taxon_rank"],function(x) which(x==taxRankPBDB))		
		#drop taxa below specified rank
		data1<-data1[rank1<=numTaxonRank,]
		#also recreate numTaxonRank
		numTaxonRank<-sapply(data1[,"taxon_rank"],function(x) which(x==taxRankPBDB))		
		#
		#create lookup table for taxon names
		taxonNameTable<-cbind(as.numeric(data1[,"taxon_no"]),as.character(data1[,"taxon_name"]))
		#add parents not listed
		parentFloat<-unique(data1[,"parent_no"])
		parentFloat<-parentFloat[is.na(match(parentFloat,taxonNameTable[,1]))]
		taxonNameTable<-rbind(taxonNameTable,cbind(parentFloat,paste("ID:",as.character(parentFloat))))
		#DONE
		#
		#now need to put together parentChild table
		#first, get table of all parentChild relationships in data
		pcAll<-cbind(as.numeric(data1[,"parent_no"]),as.numeric(data1[,"taxon_no"]))
		#then start creating final matrix: first, those of desired rank
		pcMat<-pcAll[rank1==numTaxonRank,]
		#identify IDs of parents floating without ancestors of their own
		getFloat<-function(pcDat){unique(pcDat[sapply(pcDat[,1],function(x) all(x!=pcDat[,2])),1])}
		floaters<-getFloat(pcDat=pcMat)
		#
		nCount<-0	#start tracing back
		while(length(floaters)>1){	#so only ONE root can float
			nCount<-nCount+1
			#get new relations: will 'anchor' the floaters
			anchors<-match(floaters,pcAll[,2])
			pcMat<-rbind(pcMat,pcAll[anchors[!is.na(anchors)],])	#bind to pcMat
			floatersNew<-getFloat(pcDat=pcMat)	#recalculate float
			#stopping condition, as this is a silly while() loop...
			if(length(floatersNew)>1 & identical(sort(floaters),sort(floatersNew))){
				if(!is.null(solveMissing)){
					if(solveMissing=="queryPBDB"){
						floatData<-queryMissingParents(taxaID=floatersNew,APIversion=APIversion)	
						#update taxon names in taxonNameTable
						whichUpdate<-match(floatData[,"taxon_no"],taxonNameTable[,1])
						taxonNameTable[whichUpdate[!is.na(whichUpdate)],2]<-floatData[
							!is.na(whichUpdate),"taxon_name"]
						#add any new parent taxa to taxonNameTable
						parentFloat<-unique(floatData[,"parent_no"])
						parentFloat<-parentFloat[is.na(match(parentFloat,taxonNameTable[,1]))]
						if(length(parentFloat)>0){
							taxonNameTable<-rbind(taxonNameTable,
								cbind(parentFloat,paste("ID:",as.character(parentFloat))))
							}
						#update parentChildMat, parentChildAll
						newEntries<-floatData[,c("parent_no","taxon_no")]
						pcMat<-rbind(pcMat,newEntries)
						pcAll<-rbind(pcAll,newEntries)
						}
					if(solveMissing=="mergeRoots"){
						pcMat<-rbind(pcMat,cbind("ArtificialRoot",floaters))
						taxonNameTable<-rbind(taxonNameTable,c("ArtificialRoot","ArtificialRoot"))
						message(paste0(
							"Multiple potential root-taxa artificially merged at a common root: \n",
							paste0(taxonNameTable[match(floaters,taxonNameTable[,1]),2]
								,collapse=", ")))
						}
					#regardless of method used, recalculate floaters
					floaters<-getFloat(pcDat=pcMat)
				}else{
					stop(paste0("Provided PBDB Dataset does not appear to have a \n",
						" monophyletic set of parent-child relationship pairs. \n",
						"Multiple taxa appear to be listed as parents, but are not \n",
						"listed themselves so have no parents listed: \n",
						paste0(taxonNameTable[match(floaters,taxonNameTable[,1]),2]
							,collapse=", ")))
					}
			}else{
				floaters<-floatersNew
				}
			}
		pcMat<-apply(pcMat,2,as.character)
		tree<-parentChild2taxonTree(parentChild=pcMat,tipSet=tipSet,cleanTree=cleanTree)
		#convert tip.label and node.label to taxon names from taxonNameTable
		tree$tip.label<-taxonNameTable[match(tree$tip.label,taxonNameTable[,1]),2]
		tree$node.label<-taxonNameTable[match(tree$node.label,taxonNameTable[,1]),2]
		tree$parentChild<-pcMat
		}
	#
	if(method=="Linnean"){
		#Check if show=phylo was used
		if(!any(colnames(data1)=="family")){stop("Data must be a taxonomic download with show=phylo for method='Linnean'")}
		#message that tipSet and solveMissing are ignored
		message("Linnean taxon-tree option selected, arguments 'tipSet', 'solveMissing' ignored")
		#now check and return an error if duplicate taxa of selected rank
		nDup<-sapply(nrow(data1),function(x) sum(data1[,"taxon_name"]==data1[x,"taxon_name"])>1 & data1[x,"taxon_rank"]==rank)
		if(any(nDup)){
			stop(paste0("Duplicate taxa of selected rank: ",paste0("(",which(nDup),") ",data1[nDup,"taxon_name"],collapse=", ")))}
		#filter on rank
		data1<-data1[data1[,"taxon_rank"]==rank,]
		#
		#get the fields you want
		taxonFields<-c("kingdom","phylum","class","order","family",
			"taxon_name")
		taxonData<-data1[,taxonFields]
		taxonData<-apply(taxonData,2,as.character)
		tree<-taxonTable2taxonTree(taxonTable=taxonData,cleanTree=cleanTree)
		tree$taxonTable<-taxonData
		}
	return(tree)
	}

#hidden function, don't import
translatePBDBtaxa<-function(data){
	# Do some translation
	#need to replace any empty string values with NAs (due perhaps to use of read.csv with the API)
	data[data==""]<-NA
	#if com vocab
	if(any("rnk"==colnames(data))){	
		#apparently it doesn't matter if these columns *are* present or not
		colnames(data)[colnames(data)=="acn"]<-"accepted_name"
		colnames(data)[colnames(data)=="snp"]<-"senpar_no"
		colnames(data)[colnames(data)=="rnk"]<-"taxon_rank"
		colnames(data)[colnames(data)=="nam"]<-"taxon_name"
		colnames(data)[colnames(data)=="fml"]<-"family"
		colnames(data)[colnames(data)=="odl"]<-"order"
		colnames(data)[colnames(data)=="cll"]<-"class"	
		colnames(data)[colnames(data)=="phl"]<-"phylum"	
		colnames(data)[colnames(data)=="kgl"]<-"kingdom"
		colnames(data)[colnames(data)=="par"]<-"parent_no"
		colnames(data)[colnames(data)=="oid"]<-"taxon_no"
		# taxon rank translation vectors for compact vocab
		taxRankPBDB<-c("subspecies","species","subgenus","genus","subtribe","tribe","subfamily",
			"family","superfamily","infraorder","suborder","order","superorder","infraclass",
			"subclass","class","superclass","subphylum","phylum","superphylum","subkingdom",
			"kingdom","unranked clade","informal")
		taxRankCOM<-2:26
		#change contents of "identified_rank" and "accepted_rank"
		data$taxon_rank<-sapply(data$taxon_rank,function(x) taxRankPBDB[x==taxRankCOM])
		message("compact vocab detected, relevant fields will be translated")
		}
	#if 1.1 and vocab is pbdb
	if(any(colnames(data)=="rank")){
		colnames(data)[colnames(data)=="rank"]<-"taxon_rank"
		}
	#
	#if 1.2 and there is an accepted_name column..
	if(any(colnames(data)=="accepted_name")){
		#fill empty accepted_name values with taxon_name
		nameFormal<-data[,"accepted_name"]
		nameFormal[is.na(nameFormal)]<-as.character(data[is.na(nameFormal),"taxon_name"])
		#replace taxon_name
		data[,"taxon_name"]<-nameFormal
		#
		#replace taxon_no with accepted_no
		taxNum<-data[,"accepted_no"]
		taxNum[is.na(taxNum)]<-as.character(data[is.na(taxNum),"taxon_no"])	
		data[,"taxon_no"]<-taxNum			
		#
		#also replace parent_no in the same way with senpar_no
		parNum<-data[,"senpar_no"]
		parNum[is.na(parNum)]<-as.character(data[is.na(parNum),"parent_no"])	
		data[,"parent_no"]<-parNum		
		}
	#
	return(data)
	}

#another hidden function
queryMissingParents<-function(taxaID, APIversion="1.1"){
	#drop Eukarya, as it won't return if status=senior under 1.1
	#taxaID<-as.numeric(taxaID[taxaID!="1"])
	#let's get some taxonomic data
	floatData<-read.csv(paste0("http://paleobiodb.org/",
		"data",APIversion,
		"/taxa/list.txt?id=",paste0(taxaID,collapse=","),
		"&rel=self&status=senior&vocab=pbdb"),
		stringsAsFactors=FALSE)
	if(nrow(floatData)==0){
		stop(paste("Current PBDB API would not return info for the following taxon IDs: \n",
			paste0(taxaID,collapse=", ")))}
	floatData<-translatePBDBtaxa(floatData)
	#parse down to just taxon_name, taxon_no, parent_no
	floatData<-cbind("taxon_name"=as.character(floatData[,"taxon_name"]),
		"parent_no"=as.numeric(floatData[,"parent_no"]),
		"taxon_no"=as.numeric(floatData[,"taxon_no"]))
	return(floatData)
	}
