#' Convert Simulated Taxon Data into a Cladogram
#' 
#' Convert ancestor-descendant relationships of taxa into an 'ideal' unscaled
#' cladogram, where taxa that could share true synapomorphies are arranged into
#' nested clades.
#' 
#' @details
#' This function simulates an ideal cladistic process, where the relationships
#' of a set of morphologically static taxa is resolved into a set of nested
#' hierarchial relationships (a standard cladogram), as much as would be
#' expected given the input relationships among those taxa. taxa2cladogram uses
#' information on the ancestor-descendant relationships of a bunch of taxa and
#' constructs an unscaled cladogram of the hierarcially-nesting relationships
#' among those taxa. There's no actual cladistics going on, this is just a
#' simulation of that process. If there is any chance that a set of taxa could
#' be resolved into a set of nested relationships given their
#' ancestor-descendant relationships, they will be resolved so in the output of
#' taxa2cladogram. No morphological characters are considered, we just assume
#' that if there is a nesting relationship, then it could be resolved as such.
#' This makes it the "ideal" cladogram of a simulated clade.
#' 
#' The result will probably not be fully resolved, as including both ancestor
#' and descendant taxa will generally make it impossible to produce a fully
#' nesting system of relationships. For example, consider a set of three
#' morphologically-static taxa where the first is an ancestor (either direct or
#' indirect, ala Foote, 1996) of both the second and third. If we imagine an
#' ideal cladistic analysis of the morphological characters of those three
#' taxa, this set of taxa will be unable to be broken up into
#' bifurcating-nested relationships and thus result in a polytomy. Any set of
#' ancestor-descendant relationships will have many of these, as some ancestors
#' must have more than one descendant for the clade to diversify, as noted by
#' Wagner and Erwin, 1995.
#' 
#' If there are cryptic taxa present in the output from \code{simFossilRecord}, these
#' and any of their morphologically distinguishable descendants are collapsed
#' into a polytomy to simulate the expected pattern of lack of phylogenetic
#' resolution. In addition to this merging, cryptic taxa can be dropped via the
#' argument drop.cryptic, such that only the first 'species' of each cryptic
#' taxon assemblage is listed among the tip taxa (what we would actually expect
#' to obtain, as wouldn't recognize cryptic taxa as different OTUs). By
#' default, cryptic taxa are not dropped so that the same number of taxa as in
#' the simulated data is retained.
 
#' @param taxad A five-column matrix of taxonomic data, as output by
#' \code{simFossilRecord}, after transformation with function \code{fossilRecord2fossilTaxa}.

#' @param drop.cryptic Should cryptic species be dropped (except for the
#' first)? Not dropped by default.

#' @param plot Should the result be plotted?

#' @return The resulting phylogeny without branch lengths is output as an
#' object of class phylo.
#'
#' The tip labels are the rownames from the simulation input; see documentation
#' for \code{simFossilRecord} and \code{fossilRecord2fossilTaxa} documentation for details.

#' @author David W. Bapst

#' @seealso \code{\link{simFossilRecord}}, \code{\link{taxa2phylo}}, \link{fossilRecord2fossilTaxa}

#' @references Foote, M. 1996 On the Probability of Ancestors in the Fossil
#' Record. \emph{Paleobiology} \bold{22}(2):141-151.
#' 
#' Wagner, P., and D. Erwin. 1995 Phylogenetic patterns as tests of speciation
#' models. New approaches to speciation in the fossil record. Columbia
#' University Press, New York:87-122.

#' @examples
#' 
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
#' layout(1:2)
#' cladogram<-taxa2cladogram(taxa,plot=TRUE)
#' #compare the "real" time-scaled tree of taxon last occurrences (taxa2phylo) 
#'      #to the 'ideal' cladogram
#' tree<-taxa2phylo(taxa,plot=TRUE)
#' 
#' #testing with cryptic speciation
#' recordCrypt<-simFossilRecord(p=0.1, q=0.1, prop.cryptic=0.5, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxaCrypt<-fossilRecord2fossilTaxa(recordCrypt)
#' layout(1:2)
#' parOrig<-par(no.readonly=TRUE)
#' par(mar=c(0,0,0,0))
#' cladoCrypt1<-taxa2cladogram(taxaCrypt,drop.cryptic=FALSE)
#' plot(cladoCrypt1)
#' cladoCrypt2<-taxa2cladogram(taxaCrypt,drop.cryptic=TRUE)
#' plot(cladoCrypt2)
#'
#' #reset plotting
#' par(parOrig)
#' layout(1) 
#' 
#' @export taxa2cladogram
taxa2cladogram<-function(taxad,drop.cryptic=FALSE,plot=FALSE){
	#take a taxad and turn it into an unscaled cladogram
		#do this by forming the tree as newick format first
		#essentially, this algorithm works by defining clades as only those sets of taxa which start with the FAD of the first taxon
		#This could be a 'clade' of a single OTU
		#or a clade of a bunch of taxa where one is a budding ancestor and the other are budding descendants or whatever
		#the important thing is that there's only so many nodes as there are instances where morphotaxa end/start allowing for homologies
	#require(ape)
	if(any(taxad[,6]!=taxad[,1])){
		for(i in which(taxad[,1]!=taxad[,6])){
			#reset descendants of cryptic taxa so all bud off of first cryptic species		
			taxad[taxad[,2]==taxad[i,1],2]<-taxad[i,6]
			}
		}
	if(!testParentChild(parentChild=taxad[,2:1])){stop("taxad anc-desc relationships are inconsistent")}
	tlabs<-rownames(taxad)
	desc<-lapply(taxad[,1],function(x) (taxad[taxad[,2]==x,1])[!is.na(taxad[taxad[,2]==x,1])])
	ndesc<-sapply(desc,length)
	rank<-numeric(length(ndesc))
	rank[ndesc==0]<-1
	rank[rank==0]<-NA
	while(any(is.na(rank))){
		rank<-sapply(1:length(rank),function(x) ifelse(!is.na(rank[x]),rank[x],
				1+max(rank[sapply(desc[[x]],function(y) which(y==taxad[,1]))])))}
	#okay, now all taxa are ranked by their depth from the tips
	comp<-numeric(length(ndesc))
	lab<-list()
	lab[rank==1]<-tlabs[rank==1]
	comp[rank==1]<-1
	while(any(comp==0)){
		tpot<-comp==0
		tpot2<-rank==min(rank[tpot])
		tpick<-which(tpot & tpot2)[1]
		dlab<-paste(unlist(lab[desc[[tpick]]]),",",sep="",collapse="")
		lab[[tpick]]<-paste("(",dlab,tlabs[tpick],")",sep="")
		comp[tpick]<-1
		}
	tree1<-paste(lab[[1]],";",sep="")
	tree2<-read.tree(text=tree1)
	tree2<-cleanNewPhylo(tree2)
	if(drop.cryptic & any(taxad[,6]!=taxad[,1])){
		tree2<-drop.tip(tree2,tlabs[taxad[,6]!=taxad[,1]])
		tree2<-collapse.singles(tree2)
		}
	if(plot){plot(ladderize(tree2),show.tip.label=FALSE)}
	return(tree2)
	}
