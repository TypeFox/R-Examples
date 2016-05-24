#' Simulating Extinct Clades of Monophyletic Taxa
#' 
#' This function simulates the diversification of clades composed of
#' monophyletic terminal taxa, which are distinguished in a fashion completely
#' alternative to way taxa are defined in the simulation functions
#' \code{simFossilRecord}, \code{taxa2cladogram} and \code{taxa2phylo}.
#' 
#' \code{deadTree} generates a time-scaled topology for an entirely extinct clade of a
#' specific number of tip taxa. Because the clade is extinct and assumed to
#' have gone extinct in the distant past, many details of typical birth-death
#' simulators can be ignored. If a generated clade is already conditioned upon
#' the (a) that some number of taxa was reached and (b) then the clade went
#' extinct, the topology (i.e. the distribution of branching and extinction
#' events) among the branches should be independent of the actual generating
#' rate. The frequency of nodes is a simple mathematical function of the number
#' of taxa (i.e. number of nodes is the number of taxa -1) and their placement
#' should completely random, given that we generally treat birth-death
#' processes as independent Poisson processes. Thus, in terms of generating the
#' topology, this function is nothing but a simple wrapper for the ape function
#' rtree, which randomly places splits among a set of taxa using a simple
#' algorithm (see Paradis, 2012). To match the expectation of a birth-death
#' process, new branch lengths are calculated as an exponential distribution
#' with mean 1/sumRate, where sumRate represents the sum of the branching and
#' extinction rates. Although as long as both the branching rate and extinction
#' rates are more than zero, any non-ultrametric tree is possible, only when
#' the two rates are non-zer and equal to each other will there be a high
#' chance of getting an extinct clade with many tips. Any analyses one could do
#' on a tree such as this will almost certainly give estimates of equal
#' branching and extinction rates, just because all taxa are extinct.
#' 
#' \code{simTermTaxa} produces 'terminal-taxon' datasets; datasets of clades where the
#' set of distinguishable taxa are defined as intrinsically monophyletic. (In
#' version 1.6, I referred to this as the 'candle' mode, so named from the
#' 'candling' horticultural practice and the visual conceptualization of the
#' model.) On theoretical terms, terminal-taxa datasets are what would occur if
#' (a) only descendant lineages can be sample and (b) all taxa are immediately
#' differentiated as of the last speciation event and continue to be so
#' differentiated until they go extinct. In practice, this means the taxa on
#' such a tree would represent a sample of all the terminal branches, which
#' start with some speciation event and end in an extinction event. These are
#' taken to be the true original ranges of these taxa. No further taxa can be
#' sampled than this set, whatsoever. Note that the differentiation here is a
#' result of a posteriori consideration of the phylogeny: one can't even know
#' what lineages could be sampled or the actual start points of such taxa until
#' after the entire phylogeny of a group of organisms is generated.
#' 
#' Because all evolutionary history prior to any branching events is unsampled,
#' this model is somewhat agnostic about the general model of differentiation
#' among lineages. The only thing that can be said is that synapomorphies are
#' assumed to be potentially present along every single branch, such that in an
#' ideal scenario every clade could be defined. This would suggest very high
#' anagenesis or bifurcation.
#' 
#' Because the set of observable taxa is a limited subset of the true evolution
#' history, the true taxon ranges are not a faithful reproduction of the true
#' diversity curve. See an example below.
#' 
#' \code{simTermTaxa} uses \code{deadTree} to make a phylogeny, so the only datasets produced
#' are of extinct clades. \code{simTermTaxaAdvanced} is an alternative to \code{simTermTaxa}
#' which uses \code{simFossilRecord} to generate the underlying pattern of evolutionary
#' relationships and not \code{deadTree}. The arguments are thus similar to
#' \code{simFossilRecord}, with some differences (as \code{simTermTaxaAdvanced}
#' originally called the deprecated function \code{simFossilTaxa}).
#' In particular, \code{simTermTaxaAdvanced} can be used to produce
#' simulated datasets which have extant taxa. 
#' 
#' \code{trueTermTaxaTree} is analagous to the function of \code{taxa2phylo}, in that it
#' outputs the time-scaled-phylogeny for a terminal-taxon dataset for some
#' times of observations. Unlike with the use of \code{taxa2phylo} on the output on
#' \code{simFossilRecord} (via \code{fossilRecord2fossilTaxa},
#' there is no need to use \code{trueTermTaxaTree} to obtain the true
#' phylogeny when times of extinction are the times of observation; just get
#' the \code{$tree} element from the result output by \code{simTermTaxa}.
#' 
#' Also unlike with \code{taxa2phylo}, the cladistic topology of relationships among
#' morphotaxa never changes as a function of time of observation. For obtaining
#' the 'ideal cladogram' of relationships among the terminal taxa, merely take
#' the $tree element of the output from \code{simtermTaxaData} and remove the branch
#' lengths (see below for an example).
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' @aliases termTaxa candleTaxa simCandleTaxa trueCandle simTermTaxa
#' simTermTaxaAdvanced trueTermTaxaTree deadTree

#' @param ntaxa Number of monophyletic 'terminal' taxa (tip terminals) to be
#' included on the simulated tree

#' @param sumRate The sum of the instantaneous branching and extinction rates;
#' see below.

#' @param p Instantaneous rate of speciation/branching.

#' @param q Instantaneous rate of extinction.

#' @param mintaxa Minimum number of total taxa over the entire history of a
#' clade necessary for a dataset to be accepted.

#' @param maxtaxa Maximum number of total taxa over the entire history of a
#' clade necessary for a dataset to be accepted.

#' @param mintime Minimum time units to run any given simulation before
#' stopping.

#' @param maxtime Maximum time units to run any given simulation before
#' stopping.

#' @param minExtant Minimum number of living taxa allowed at end of
#' simulations.

#' @param maxExtant Maximum number of living taxa allowed at end of
#' simulations.

#' @param min.cond If TRUE, the default, simulations are stopped when they meet
#' all minimum conditions. If FALSE, simulations will continue until they hit
#' maximum conditions, but are only accepted as long as they still meet all
#' minimum conditions in addition.

#' @param TermTaxaRes The list output produced by simTermTaxa

#' @param time.obs A per-taxon vector of times of observation for the taxa in
#' TermTaxaRes

#' @return \code{deadTree} gives time-scaled phylo object, with a $root.time element.
#' As discussed above, the result is always an extinct phylogeny of exactly
#' \code{ntaxa}.
#' 
#' \code{simTermTaxa} and \code{simTermTaxaAdvanced} both produce a list with two components:
#' \code{$taxonRanges} which is a two-column matrix where each row gives the true
#' first and last appearance of observable taxa and \code{$tree} which is a
#' time-scaled phylogeny with end-points at the true last appearance time of
#' taxa.
#' 
#' \code{trueTermTaxaTree} produces a time-scaled tree as a phylo object, which
#' describes the relationships of populations at the times of observation given
#' in the time.obs argument.

#' @author David W. Bapst

#' @seealso deadtree is simply a wraper of the function \code{rtree} in ape.
#' 
#' For a very different way of simulating diversification in the fossil record,
#' see \code{\link{simFossilRecord}}, \code{\link{fossilRecord2fossilTaxa}},
#' \code{\link{taxa2phylo}} and \code{\link{taxa2cladogram}}.

#' @references Paradis, E. (2012) \emph{Analysis of Phylogenetics and Evolution
#' with R (Second Edition).} New York: Springer.

#' @examples
#' 
#' set.seed(444)
#' #example for 20 taxa
#' termTaxaRes<-simTermTaxa(20)
#' 
#' #let look at the taxa...
#' taxa<-termTaxaRes$taxonRanges
#' taxicDivCont(taxa)
#' #because ancestors don't even exist as taxa
#' 	#the true diversity curve can go to zero
#' 	#kinda bizarre!
#' 
#' #the tree should give a better idea
#' tree<-termTaxaRes$tree
#' phyloDiv(tree)
#' #well, okay, its a tree. 
#' 
#' #get the 'ideal cladogram' ala taxa2cladogram
#'     #much easier with terminal-taxa simulations as no paraphyletic taxa
#' cladogram<-tree
#' cladogram$edge.length<-NULL
#' plot(cladogram)
#' 
#' #trying out trueTermTaxaTree
#' #random times of observation: uniform distribution
#' time.obs<-apply(taxa,1,function(x) runif(1,x[2],x[1]))
#' tree1<-trueTermTaxaTree(termTaxaRes,time.obs)
#' layout(1:2)
#' plot(tree)
#' plot(tree1)
#' 
#' layout(1)
#' 
#' #let's look at the change in the terminal branches
#' plot(tree$edge.length,tree1$edge.length)
#' #can see some edges are shorter on the new tree, cool
#' 
#' #let's now simulate sampling and use FADs
#' layout(1:2)
#' plot(tree);axisPhylo()
#' FADs<-sampleRanges(termTaxaRes$taxonRanges,r=0.1)[,1]
#' tree1<-trueTermTaxaTree(termTaxaRes,FADs)
#' plot(tree1);axisPhylo()
#' 
#' #can condition on sampling some average number of taxa
#' #analagous to deprecated function simFossilTaxa_SRcond
#' r<-0.1
#' avgtaxa<-50
#' sumRate<-0.2
#' #avg number necc for an avg number sampled
#' ntaxa_orig<-avgtaxa/(r/(r+sumRate))	
#' termTaxaRes<-simTermTaxa(ntaxa=ntaxa_orig,sumRate=sumRate)
#' #note that conditioning must be conducted using full sumRate
#' #this is because durations are functions of both rates
#' #just like in bifurcation
#' 
#' #use advanced version of simTermTaxa: simTermTaxaAdvanced
#'     #allows for extant taxa in a term-taxa simulation
#' 
#' #with min.cond
#' termTaxaRes<-simTermTaxaAdvanced(p=0.1,q=0.1,mintaxa=50,
#'     maxtaxa=100,maxtime=100,minExtant=10,maxExtant=20,min.cond=TRUE)
#' #notice that arguments are similar to simFossilRecord
#' 	# and somewhat more similar to deprecated function simFossilTaxa ;P
#' plot(termTaxaRes$tree)
#' Ntip(termTaxaRes$tree)
#' 
#' #without min.cond
#' termTaxaRes<-simTermTaxaAdvanced(p=0.1,q=0.1,mintaxa=50,
#'     maxtaxa=100,maxtime=100,minExtant=10,maxExtant=20,min.cond=FALSE)
#' plot(termTaxaRes$tree)
#' Ntip(termTaxaRes$tree)
#' 
#' layout(1)
#' 
#' @name termTaxa
#' @rdname termTaxa
#' @export
simTermTaxa<-function(ntaxa,sumRate=0.2){
		#previously known as "candle"
	#This function will make ideal cladist datasets
	#all taxa will be descendants
		#all evolution prior to branching or differentiation is unsampled
		#taxa do not differentiatiate over those times
	#require(ape)
	tree<-deadTree(ntaxa=ntaxa,sumRate=sumRate)
	#taxonNames<-tree$tip.label
	termEdge<-sapply(tree$edge[,2],function(x) any(x==(1:ntaxa)))
	#termAnc<-tree$edge[termEdge,1]
	taxonDurations<-tree$edge.length[termEdge]
	nodeDist<-node.depth.edgelength(tree)
	taxonLADs<-tree$root.time-nodeDist[1:ntaxa]
	taxonFADs<-taxonLADs+taxonDurations
	taxonRanges<-cbind(taxonFADs,taxonLADs)
	rownames(taxonRanges)<-tree$tip.label[tree$edge[termEdge,2]]
	res<-list(taxonRanges=taxonRanges,tree=tree)
	return(res)
	}

#' @rdname termTaxa
#' @export
simTermTaxaAdvanced<-function(p=0.1,q=0.1,mintaxa=1,maxtaxa=1000,mintime=1,maxtime=1000,
		minExtant=0,maxExtant=NULL,min.cond=TRUE){
		#previously known as "candle"
	#This function will make ideal cladist datasets
	#all taxa will be monophyletic descendants
		#all evolution prior to branching or differentiation is unsampled
		#taxa do not differentiate over those times
	#extant example
		#p=0.1;q=0.1;mintaxa=50;maxtaxa=100;mintime=1;maxtime=1000;minExtant=10;maxExtant=20;min.cond=FALSE
	#require(ape)
	#taxa<-simFossilTaxa(p=p,q=q,mintaxa=mintaxa,maxtaxa=maxtaxa,mintime=mintime,maxtime=maxtime,
	#	minExtant=minExtant,maxExtant=maxExtant,min.cond=min.cond,nruns=1,
	#	anag.rate=0,prop.bifurc=0,prop.cryptic=0,count.cryptic=FALSE,print.runs=FALSE,plot=FALSE)
	record<-simFossilRecord(p=p, q=q, r = 0, nruns = 1,
		nTotalTaxa=c(mintaxa,maxtaxa),totalTime=c(mintime,maxtime),nExtant=c(minExtant,maxExtant),
		anag.rate = 0, prop.bifurc = 0, prop.cryptic = 0, modern.samp.prob = 1, startTaxa = 1, 
		nSamp = c(0, 1000), tolerance = 10^-4, maxStepTime = 0.01,
		shiftRoot4TimeSlice = "withExtantOnly", count.cryptic = FALSE,
		negRatesAsZero = TRUE, print.runs = FALSE, sortNames = FALSE,
		plot = FALSE)
	taxa<-fossilRecord2fossilTaxa(record)
	tree<-dropZLB(taxa2phylo(taxa))
	ntaxa<-Ntip(tree)
	#taxonNames<-tree$tip.label
	termEdge<-sapply(tree$edge[,2],function(x) any(x==(1:ntaxa)))
	termNodes<-tree$edge[termEdge,2]
	#termAnc<-tree$edge[termEdge,1]
	taxonDurations<-tree$edge.length[termEdge]
	nodeDist<-node.depth.edgelength(tree)
	taxonLADs<-tree$root.time-nodeDist[termNodes]
	taxonFADs<-taxonLADs+taxonDurations
	taxonRanges<-cbind(taxonFADs,taxonLADs)
	rownames(taxonRanges)<-tree$tip.label[termNodes]
	res<-list(taxonRanges=taxonRanges,tree=tree)
	return(res)
	}

#' @rdname termTaxa
#' @export
trueTermTaxaTree<-function(TermTaxaRes,time.obs){
		#for a term-taxa (candle) tree datasets
	#given a per-taxon vector of observation times, returns the true tree
	#time.obs must have taxon names
	#require(ape)
	#first, check if the times of observations are outside of original taxon ranges
	taxR<-TermTaxaRes$taxonRanges
	nameMatch<-match(names(time.obs),rownames(taxR))
	if(any(is.na(nameMatch))){stop("ERROR: names on time.obs and in TermTaxaRes don't match")}
	if(any(sapply(1:length(time.obs),function(x)
			if(is.na(time.obs[x])){
				FALSE
			}else{
				(time.obs[x]>taxR[nameMatch[x],1])|(time.obs[x]<taxR[nameMatch[x],2])
				}))
		){
			stop("ERROR: Given time.obs are outside of the original taxon ranges")}
	#now onwards with the actual function
	tree1<-TermTaxaRes$tree
	#termEdge<-sapply(tree1$edge[,2],function(x)
	#	any(x==(1:Ntip(tree1))))	
	newDurations<-taxR[nameMatch,1]-time.obs
	if(is.null(names(time.obs))){
		stop("ERROR: No taxon names on observation vector?")
		}
	tipMatch<-sapply(1:Ntip(tree1),function(x)
		which(tree1$tip.label[x]==names(time.obs)))
	dropTaxa<-character()
	for(i in 1:Ntip(tree1)){
		newDur<-newDurations[tipMatch[i]]
		if(!is.na(newDur)){
			if(tree1$edge.length[tree1$edge[,2]==i]<newDur){
				stop("New duration longer than original taxon ranges?")
				}
			tree1$edge.length[tree1$edge[,2]==i]<-newDur
		}else{
			dropTaxa<-c(dropTaxa,tree1$tip.label[i])
			}
		}
	treeNoDrop<-tree1		#the root.time can't change cause the term branches got shifted
	if(length(dropTaxa)>0){
		tree1<-drop.tip(tree1,dropTaxa)
		}
	#need to correct root.time if basal outgroups were removed
	tree1<-fixRootTime(treeNoDrop,tree1)	#use the tree after adjusting the term branch lengths
	return(tree1)
	}

#' @rdname termTaxa
#' @export
deadTree<-function(ntaxa,sumRate=0.2){
	#all taxa will always be extinct
		#assuming the clade is extinct, p and q are meaningless
		#all trees would have us assume that p=q if analyzed...
		#we'll just have a 'rate' of events (branching or ext; p+q)
		#the probability of a branch ending in branching or ext is meaningless!
			#by definition equal proportion must branch and go extinct
			#anything else should be indep in a b-d process
	#lets make the phylogeny 
	#require(ape)
	tree<-rtree(ntaxa)
	tree$edge.length<-rexp(ntaxa+ntaxa-2,sumRate)
	tree$root.time<-max(node.depth.edgelength(tree))+200
	return(tree)
	}
