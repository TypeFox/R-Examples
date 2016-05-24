#' Simulate Cladogenetic Trait Evolution
#' 
#' This function simulates trait evolution at each speciation/branching event
#' in a matrix output from \code{simFossilRecord}, after transformation with
#' \code{fossilRecord2fossilTaxa}.
#' 
#' @details This function simulates continuous trait evolution where change occurs under
#' a Brownian model, but only at events that create new distinct morphotaxa
#' (i.e. species as recognized in the fossil record), either branching events
#' or anagenesis (pseudospeciation). These are the types of morphological
#' differentiation which can be simulated in the function \code{simFossilRecord}. This
#' is sometimes referred to as cladogenetic or speciation trait evolution and
#' is related to Puncuated Equilibrium theory. Anagenestic shifts aren't
#' cladogenetic events per se (no branching!), so perhaps the best way to this
#' of this function is it allows traits to change anytime \code{simFossilRecord} created
#' a new 'morphotaxon' in a simulation.
#' 
#' Importantly, trait changes only occur at the base of 'new' species, thus
#' allowing cladogenetic trait evolution to be asymmetrical at branching
#' points: i.e. only one branch actually changes position in trait-space, as
#' expected under a budding cladogenesis model. This distinction is important
#' as converting the taxa matrix to a phylogeny and simulating the trait
#' changes under a 'speciational' tree-transformation would assume that
#' divergence occurred on both daughter lineages at each node. (This has been
#' the standard approach for simulating cladogenetic trait change on trees).
#' 
#' Cryptic taxa generated with \code{prop.cryptic} in \code{simFossilRecord} will not differ at
#' all in trait values. These species will all be identical.
#' 
#' See this link for additional details:
#' http://nemagraptus.blogspot.com/2012/03/simulating-budding-cladogenetictrait.html
#' 
#' @param taxa A five-column matrix of taxonomic data, as output by
#' \code{fossilRecord2fossilTaxa} after simulation with \code{simFossilRecord}
#' (or via the deprecated function \code{simFossilTaxa})

#' @param rate rate of trait change; variance of evolutionary change
#' distribution per speciation event

#' @param meanChange Mean change per speciation event. Default is 0; change to
#' simulate 'active' speciational trends, where the expected change at each
#' speciational event is non-zero.

#' @param rootTrait The trait value of the first taxon in the dataset; set to 0
#' by default.

#' @return Retuns a vector of trait values for each taxon, with value names
#' being the taxa IDs (column 1 of the input) with a 't' pasted (as with rtree
#' in the ape library).

#' @author David W. Bapst

#' @seealso \code{\link{simFossilRecord}},
#' 
#' This function is similar to Brownian motion simulation functions such as
#' \code{rTraitCont} in ape, \code{sim.char} in geiger and \code{fastBM} in
#' phytools.
#' 
#' See also \code{\link{unitLengthTree}} in this package and
#' \code{speciationalTree} in the package geiger. These are tree transformation
#' functions; together with BM simulation functions, they would be expected to
#' have a similar effect as this function (when cladogenesis is 'bifurcating'
#' and not 'budding'; see above).

#' @examples
#' 
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,1000), plot=TRUE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' trait <- cladogeneticTraitCont(taxa)
#' tree <- taxa2phylo(taxa)
#' plotTraitgram(trait,tree,conf.int=FALSE)
#' 
#' #with cryptic speciation
#' record<-simFossilRecord(p=0.1, q=0.1, prop.cryptic=0.5, 
#'	nruns=1, nTotalTaxa=c(30,1000), plot=TRUE)
#' taxa<-fossilRecord2fossilTaxa(record)
#' trait <- cladogeneticTraitCont(taxa)
#' tree <- taxa2phylo(taxa)
#' plotTraitgram(trait,tree,conf.int=FALSE)
#' 

#' @export cladogeneticTraitCont
cladogeneticTraitCont<-function(taxa,rate=1,meanChange=0,rootTrait=0){
	#simulate speciational trait evolution for datasets from simFossilRecord
	#idiot proofing
	taxa<-taxa[order(taxa[,1]),]
	if(any(taxa[-1,2]>=taxa[-1,1])){stop("Ancestors have higher IDs than Descendants?")}
	anctaxa<-sapply(taxa[-1,2],function(x) which(x==taxa[,1]))
	traits<-rootTrait
	for(i in 2:nrow(taxa)){
		if(taxa[i,1]==taxa[i,6]){
			traits[i]<-traits[anctaxa[i-1]]+rnorm(1,mean=meanChange,sd=sqrt(rate))
		}else{
			traits[i]<-traits[anctaxa[i-1]]
			}
		}
	names(traits)<-rownames(taxa)
	return(traits)
	}
