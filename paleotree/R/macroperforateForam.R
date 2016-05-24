#' Ancestor-Descendant Relationships for Macroperforate Foraminifera, from Aze et al. (2011)
#'
#' An example dataset of ancestor-descendent relationships and first and last appearance dates for
#' a set of macroperforate Foramanifera, taken from the supplemental materials of Aze et al. (2011).
#' This dataset is included here primarily for testing functions \code{parentChild2taxonTree}
#' and \code{taxa2phylo}.

#' @name macroperforateForam

#' @rdname macroperforateForam

#' @aliases foramAL foramAM foramALb foramAMb

#' @details
#' This example dataset is composed of four tables, each containing information
#' on the ancestor-descendant relationships and first and last appearances of
#' species of macroperforate foraminifera species from the fossil record.
#' Each of the four tables are for the same set of taxa, but divide and
#' concatanate the included foram species in four different ways, relating to
#' the use of morpospecies versus combined anagenetic lineages (see Ezard et
#' al., 2012), and whether taxa are retained as units related by budding-cladogensis
#' or the splitting of taxa at branching points to create a fully 'bifurcating' set
#' of relationships, independent of ancestral morphotaxon persistance through branching
#' events. See the examples section for more details.

#' @format 
#' The 'foramAM' and 'foramAL' tables include budding taxon units
#' for morphospecies and lineages respective, with four columns:
#' taxon name, ancestral taxon's name, first appearance date and last appearance
#' date (note that column headings vary). The 'foramAMb' and 'foramALb' tables are
#' composed of data for the same taxon units as the previous
#" set, except parent taxa that persist through
#' branching events are split so that the relationships are fully 'bifurcating', rather
#' than 'budding'. As this obscures taxonomic identity, taxon identification labels
#' are included in an additional, fifth column in these tables.
#' See the examples section for more details.

#' @source 
#' This dataset is obtained from the supplementary materials of, specifically
#' 'Appendix S5':
#'
#' Aze, T., T. H. G. Ezard, A. Purvis, H. K. Coxall, D. R. M. Stewart,
#'  B. S. Wade, and P. N. Pearson. 2011. A phylogeny of Cenozoic macroperforate
#' planktonic foraminifera from fossil data. \emph{Biological Reviews} 86(4):900-927.

#' @references
#'
#' This dataset has been used or referenced in a number of works, including:
#'
#' Aze, T., T. H. G. Ezard, A. Purvis, H. K. Coxall, D. R. M. Stewart, B. S. Wade, and P. N. Pearson. 2013.
#' Identifying anagenesis and cladogenesis in the fossil record.
#' \emph{Proceedings of the National Academy of Sciences} 110(32):E2946-E2946.
#' 
#' Ezard, T. H. G., T. Aze, P. N. Pearson, and A. Purvis. 2011. Interplay Between Changing Climate and Species'
#' Ecology Drives Macroevolutionary Dynamics. \emph{Science} 332(6027):349-351.
#' 
#' Ezard, T. H. G., P. N. Pearson, T. Aze, and A. Purvis. 2012. The meaning of birth and death (in
#' macroevolutionary birth-death models). \emph{Biology Letters} 8(1):139-142.
#' 
#' Ezard, T. H. G., G. H. Thomas, and A. Purvis. 2013. Inclusion of a near-complete fossil record reveals
#' speciation-related molecular evolution. \emph{Methods in Ecology and Evolution} 4(8):745-753.
#' 
#' Strotz, L. C., and A. P. Allen. 2013. Assessing the role of cladogenesis in macroevolution by integrating
#' fossil and molecular evidence. \emph{Proceedings of the National Academy of Sciences} 110(8):2904-2909.
#' 
#' Strotz, L. C., and A. P. Allen. 2013. Reply to Aze et al.: Distinguishing speciation modes based on
#' multiple lines of evidence. \emph{Proceedings of the National Academy of Sciences} 110(32):E2947-E2947.
#' 

#' @keywords datasets

#' @docType data

#' @examples
#' 
#' # Following Text Reproduced from Aze et al. 2011's Supplemental Material
#' # Appendix S5
#' # 
#' # 'Data required to produce all of the phylogenies included in the manuscript
#' # using paleoPhylo (Ezard & Purvis, 2009) a free software package to draw
#' # paleobiological phylogenies in R.'
#' #
#' # 'The four tabs hold different versions of our phylogeny:
#' #	 aMb: fully bifurcating morphospecies phylogeny
#' #	 aM: budding/bifurcating morphospecies phylogeny
#' #	 aLb: fully bifurcating lineage phylogeny
#' #	 aL: budding/bifurcating lineage phylogeny
#' #
#' # 'Start Date gives the first occurence of the species according
#' # to the particular phylogeny; End Date gives the last occurence
#' # according to the particular phylogeny.'
#' 
#' \dontrun{
#' 
#' #load the data 
#' 	#given in supplemental as XLS sheets
#' 	#converted to separate tab-deliminated text files
#' 
#' # aM: budding/bifurcating morphospecies phylogeny
#' foramAM<-read.table(file.choose(),stringsAsFactors=FALSE,header=TRUE)
#' # aL: budding/bifurcating lineage phylogeny
#' foramAL<-read.table(file.choose(),stringsAsFactors=FALSE,header=TRUE)
#' # aMb: fully bifurcating morphospecies phylogeny
#' foramAMb<-read.table(file.choose(),stringsAsFactors=FALSE,header=TRUE)
#' # aLb: fully bifurcating lineage phylogeny
#' foramALb<-read.table(file.choose(),stringsAsFactors=FALSE,header=TRUE)
#' 
#' save.image("macroperforateForam.rdata")
#' 
#' }
#' 
#' #instead, we'll just load the data directly
#' data(macroperforateForam)
#' 
#' #Two distinctions among the four datasets:
#' #(1): morphospecies vs morphospecies combined into sequences of anagenetic
#' 	# morpospecies referred to as 'lineages'. Thus far more morphospecies
#' 	# than lineages. The names of lineages are given as the sequence of
#' 	# their respective component morphospecies.
#' #(2): Datasets where taxon units (morphospecies or lineages) are broken up
#' 	# at 'budding' branching events (where the ancestral taxon persists)
#' 	# so that final dataset is 'fully bifurcating', presumably
#' 	# to make comparison easier to extant-taxon only datasets.
#' 	# (This isn't a limitation for paleotree, though!).
#' 	# This division of taxon units requires abstracting the taxon IDs,
#' 	# requiring another column for Species Name.
#' 
#' dim(foramAM)
#' dim(foramAL)
#' dim(foramAMb)
#' dim(foramALb)
#' 
#' #Need to convert these to same format as fossilRecord2fossilTaxa output.
#' 	#those 'taxa' tables has 6 columns:
#' 	#taxon.id ancestor.id orig.time ext.time still.alive looks.like
#' 
#' #for the purposes of this, we'll make taxon.id=looks.like
#' 	# (That's only for simulating cryptic speciation anyway)
#' #still.alive should be TRUE (1) if ext.time=0
#' 
#' #a function to convert Aze et al's suppmat to paleotree-readable format
#' 
#' createTaxaData<-function(table){
#' 	#reorder table by first appearance time
#' 	table<-table[order(-as.numeric(table[,3])),]
#' 	ID<-1:nrow(table)
#' 	anc<-sapply(table[,2],function(x)
#' 		if(!is.na(x)){
#' 			which(x==table[,1])
#' 		}else{ NA })
#' 	stillAlive<-as.numeric(table[,4]==0)
#' 	ages<-cbind(as.numeric(table[,3]),as.numeric(table[,4]))
#' 	res<-cbind(ID,anc,ages,stillAlive,ID)
#' 	colnames(res)<-c('taxon.id','ancestor.id','orig.time',
#' 		'ext.time','still.alive','looks.like')
#' 	rownames(res)<-table[,1]
#' 	return(res)
#' 	}
#' 
#' taxaAM<-createTaxaData(foramAM)
#' taxaAMb<-createTaxaData(foramAMb)
#' taxaAL<-createTaxaData(foramAL)
#' taxaALb<-createTaxaData(foramALb)
#' 
#' ##################################
#' 
#' #Checking Ancestor-Descendant Relationships for Irregularities
#' 
#' #For each of these, there should only be a single taxon
#' 	# without a parent listed (essentially, the root ancestor)
#' 
#' countParentsWithoutMatch<-function(table){
#'     	parentMatch<-match(unique(table[,2]),table[,1])
#'     	sum(is.na(parentMatch))
#' 	}
#' 
#' #test this on the provided ancestor-descendant relationships
#' countParentsWithoutMatch(foramAM)
#' countParentsWithoutMatch(foramAL)
#' countParentsWithoutMatch(foramAMb)
#' countParentsWithoutMatch(foramALb)
#' 
#' #and on the converted datasets
#' countParentsWithoutMatch(taxaAM)
#' countParentsWithoutMatch(taxaAL)
#' countParentsWithoutMatch(taxaAMb)
#' countParentsWithoutMatch(taxaALb)
#' 
#' \donttest{ 
#' 
#' #can construct the parentChild2taxonTree
#' 	#using the ancestor-descendant relationships 
#' 
#' #can be very slow...
#' 
#' treeAM<-parentChild2taxonTree(foramAM[,2:1])
#' treeAL<-parentChild2taxonTree(foramAL[,2:1])
#' treeAMb<-parentChild2taxonTree(foramAMb[,2:1])
#' treeALb<-parentChild2taxonTree(foramALb[,2:1])
#' 
#' layout(matrix(1:4,2,2))
#' plot(treeAM,main='treeAM',show.tip.label=FALSE)
#' plot(treeAL,main='treeAL',show.tip.label=FALSE)
#' plot(treeAMb,main='treeAMb',show.tip.label=FALSE)
#' plot(treeALb,main='treeALb',show.tip.label=FALSE)
#' 
#' # FYI 
#' # in case you were wondering
#' # you would *not* time-scale these Frankenstein monsters
#' 
#' }
#' 
#' ###########################################
#' 
#' # Checking stratigraphic ranges
#' 
#' # do all first occurrence dates occur before last occurrence dates?
#' 	# we'll check the original datasets here
#' 
#' checkFoLo<-function(data){
#' 	diffDate<-data[,3]-data[,4]	#subtract LO from FO
#' 	isGood<-all(diffDate>=0)	#is it good
#' 	return(isGood)
#' 	}
#' 
#' checkFoLo(foramAM)
#' checkFoLo(foramAL)
#' checkFoLo(foramAMb)
#' checkFoLo(foramALb)
#' 
#' #cool, but do all ancestors appear before their descendents?
#' 	# easier to check unified fossilRecord2fossilTaxa format here
#' 
#' checkAncOrder<-function(taxa){
#' 	#get ancestor's first occurrence
#' 	ancFO<-taxa[taxa[,2],3]
#' 	#get descendant's first occurrence	
#' 	descFO<-taxa[,3]
#' 	diffDate<-ancFO-descFO	#subtract descFO from ancFO
#' 	#remove NAs due to root taxon
#' 	diffDate<-diffDate[!is.na(diffDate)]
#' 	isGood<-all(diffDate>=0)	#is it all good	
#' 	return(isGood)
#' 	}
#' 
#' checkAncOrder(taxaAM)
#' checkAncOrder(taxaAL)
#' checkAncOrder(taxaAMb)
#' checkAncOrder(taxaALb)
#' 
#' #now, are there gaps between the last occurrence of ancestors
#' 	# and the first occurrence of descendents?
#' 	# (shall we call these 'stratophenetic ghost branches'?!)
#' 	# These shouldn't be problematic, but do they occur in this data?
#' # After all, fossilRecord2fossilTaxa output tables are designed for
#' 	   # fully observed simulated fossil records with no gaps.
#' 
#' sumAncDescGap<-function(taxa){
#' 	#get ancestor's last occurrence
#' 	ancLO<-taxa[taxa[,2],4]
#' 	#get descendant's first occurrence	
#' 	descFO<-taxa[,3]
#' 	diffDate<-ancLO-descFO	#subtract descFO from ancFO
#' 	#remove NAs due to root taxon
#' 	diffDate<-diffDate[!is.na(diffDate)]
#' 	#should be negative or zero, positive values are gaps
#' 	gaps<-c(0,diffDate[diffDate>0])
#' 	sumGap<-sum(gaps)
#' 	return(sumGap)
#' 	}
#' 
#' #get the total gap between ancestor LO and child FO
#' sumAncDescGap(taxaAM)
#' sumAncDescGap(taxaAL)
#' sumAncDescGap(taxaAMb)
#' sumAncDescGap(taxaALb)
#' 
#' #It appears there is *no* gaps between ancestors and their descendants
#' 	#in the Aze et al. foram dataset... wow!
#' 
#' ###############
#' 
#' \donttest{ 
#' 
#' # Creating time-scaled phylogenies from the Aze et al. data
#' 
#' # Aze et al. (2011) defines anagenesis such that taxa may overlap
#' # in time during a transitional period (see Ezard et al. 2012
#' # for discussion of this definition). Thus, we would expect that
#' # paleotree obtains very different trees for morphospecies versus
#' # lineages, but very similar phylogenies for datasets where budding
#' # taxa are retained or arbitrarily broken into bifurcating units.
#' 
#' # We can use the function taxa2phylo to directly create
#' # time-scaled phylogenies from the Aze et al. stratophenetic data
#' 
#' timetreeAM<-taxa2phylo(taxaAM)
#' timetreeAL<-taxa2phylo(taxaAL)
#' timetreeAMb<-taxa2phylo(taxaAMb)
#' timetreeALb<-taxa2phylo(taxaALb)
#' 
#' layout(matrix(1:4,2,2))
#' plot(timetreeAM,main='timetreeAM',show.tip.label=FALSE)
#' axisPhylo()
#' plot(timetreeAL,main='timetreeAL',show.tip.label=FALSE)
#' axisPhylo()
#' plot(timetreeAMb,main='timetreeAMb',show.tip.label=FALSE)
#' axisPhylo()
#' plot(timetreeALb,main='timetreeALb',show.tip.label=FALSE)
#' axisPhylo()
#' 
#' #visually compare the two pairs we expect to be close to identical
#' 
#' #morpospecies
#' layout(1:2)
#' plot(timetreeAM,main='timetreeAM',show.tip.label=FALSE)
#' axisPhylo()
#' plot(timetreeAMb,main='timetreeAMb',show.tip.label=FALSE)
#' axisPhylo()
#' 
#' #lineages
#' layout(1:2)
#' plot(timetreeAL,main='timetreeAL',show.tip.label=FALSE)
#' axisPhylo()
#' plot(timetreeALb,main='timetreeALb',show.tip.label=FALSE)
#' axisPhylo()
#' 
#' layout(1)
#' 
#' #compare the summary statistics of the trees
#' Ntip(timetreeAM)
#' Ntip(timetreeAMb)
#' Ntip(timetreeAL)
#' Ntip(timetreeALb)
#' # very different!
#' 
#' # after dropping anagenetic zero-length-terminal-edge ancestors
#' # we would expect morphospecies and lineage phylogenies to be very similar
#' 
#' #morphospecies
#' Ntip(dropZLB(timetreeAM))
#' Ntip(dropZLB(timetreeAMb))
#' #identical!
#' 
#' #lineages
#' Ntip(dropZLB(timetreeAL))
#' Ntip(dropZLB(timetreeALb))
#' # ah, very close, off by a single tip
#' # ...probably a very short ZLB outside tolerance
#' 
#' #we can create some diversity plots to compare
#' 
#' multiDiv(data=list(timetreeAM,timetreeAMb),
#' 	plotMultCurves=TRUE)
#' 
#' multiDiv(data=list(timetreeAL,timetreeALb),
#' 	plotMultCurves=TRUE)
#' 
#' # we can see that the morphospecies datasets are identical
#' 	# that's why we can only see one line
#' # some very slight disagreement between the lineage datasets
#' 	# around ~30-20 Ma
#' 
#' #can also compare morphospecies and lineages diversity curves
#' 
#' multiDiv(data=list(timetreeAM,timetreeAL),
#' 	plotMultCurves=TRUE)
#' 
#' #they are similar, but some peaks are missing from lineages
#' 	# particularly around ~20-10 Ma
#' 
#' 
#' }
#' 
#' 
#'
NULL

