#' Three Rate Calibrated 'a posteriori' Time-Scaling of Paleo-Phylogenies
#' 
#' Time-scales an unscaled cladogram of fossil taxa, using information on their
#' ranges and estimates of the instantaneous rates of branching, extinction and
#' sampling. The output is a sample of a posteriori time-scaled trees, as resulting from a
#' stochastic algorithm which samples observed gaps in the fossil record with
#' weights calculated based on the input rate estimates. This function also
#' uses the three-rate calibrated time-scaling algorithm to stochastically
#' resolve polytomies and infer potential ancestor-descendant relationships,
#' simultaneous with the time-scaling treatment.
#' 
#' @details The three-rate calibrated ("cal3") algorithm time-scales trees a posteriori by
#' stochastically picking node divergence times relative to a probability
#' distribution of expected waiting times between speciation and first
#' appearance in the fossil record. This algorithm is extended to apply to
#' resolving polytomies and designating possible ancestor-descendant
#' relationships. The full details of this method, its sister method the
#' sampling-rate-calibrated time-scaling method and the details of the algorithm
#' will be given in Bapst (in prep for Paleobiology).
#' 
#' Briefly, cal3 time-scaling is done by examining each node separately, moving
#' from the root upwards. Ages of branching nodes are constrained below by the
#' ages of the nodes below them (except the root; hence the need for the
#' root.max argument) and constrained above by the first appearance dates
#' (FADs) of the daughter lineages. The position of the branching event within
#' this constrained range implies different amounts of unobserved evolutionary
#' history. cal3 considers a large number of potential positions for the
#' branching node (the notation in the code uses the analogy of viewing the
#' branching event as a 'zipper') and calculates the summed unobserved
#' evolutionary history implied by each branching time. The probability density
#' of each position is then calculated under a gamma distribution with a shape
#' parameter of 2 (implying that it is roughly the sum of two normal waiting
#' times under an exponential) and a rate parameter which takes into account
#' both the probability of not observing a lineage of a certain duration and
#' the 'twiginess' of the branch, i.e. the probability of having short-lived
#' descendants which went extinct and never were sampled (ala Friedman and
#' Brazeau, 2011). These densities calculated under the gamma distribution are
#' then used as weights to stochastically sample the possible positions for the
#' branching node. This basic framework is extended to polytomies by allowing a
#' branching event to fall across multiple potential lineages, adding each
#' lineage one by one, from earliest appearing to latest appearing (the code
#' notation refers to this as a 'parallel zipper').
#' 
#' As with many functions in the paleotree library, absolute time is always
#' decreasing, i.e. the present day is zero.
#' 
#' These functions will intuitively drop taxa from the tree with NA for range
#' or that are missing from timeData.
#' 
#' The sampling rate used by cal3 methods is the instantaneous sampling rate,
#' as estimated by various other function in the paleotree package. See
#' \code{\link{make_durationFreqCont}} for more details.
#' If you have the per-time unit sampling
#' probability ('R' as opposed to 'r') look at the sampling parameter
#' conversion functions also included in this package
#' (e.g. \code{\link{sProb2sRate}}). Most datasets will probably use
#' \code{\link{make_durationFreqDisc}} and \code{sProb2sRate}
#' prior to using this function, as shown in an example below.
#' 
#' The branching and extinction rate are the 'per-capita' instantaneous
#' origination/extinction rates for the taxic level of the tips of the tree
#' being time-scaled. Any user of the cal3 time-scaling method has multiple
#' options for estimating these rates. One is to separately calculate the
#' per-capita rates (following the equations in Foote, 2001) across multiple
#' intervals and take the mean for each rate. A second, less preferred option,
#' would be to use the extinction rate calculated from the sampling rate above
#' (under ideal conditions, this should be very close to the mean 'per-capita'
#' rate calculated from by-interval FADs and LADs). The branching rate in this
#' case could be assumed to be very close to the extinction rate, given the
#' tight relationship observed in general between these two (Stanley, 1976; see
#' Foote et al., 1999, for a defense of this approach), and thus the extinction
#' rate estimate could be used also for the branching rate estimate. (This is
#' what is done for the examples below.) A third option for calculating all
#' three rates simultaneously would be to apply likelihood methods developed by
#' Foote (2002) to forward and reverse survivorship curves. Note that only one
#' of these three suggested methods is implemented in \code{paleotree}: estimating the
#' sampling and extinction rates from the distribution of taxon durations via
#' \code{make_durationFreqCont} and \code{make_durationFreqDisc}.
#' 
#' By default, the cal3 functions will consider that ancestor-descendant
#' relationships may exist among the given taxa, under a budding cladogenetic
#' or anagenetic modes. Which tips are designated as which is given by two
#' additional elements added to the output tree, $budd.tips (taxa designated as
#' ancestors via budding cladogenesis) and $anag.tips (taxa designated as
#' ancestors via anagenesis). This can be turned off by setting \code{anc.wt = 0}. As
#' this function may infer anagenetic relationships during time-scaling, this
#' can create zero-length terminal branches in the output. Use
#' \code{\link{dropZLB}} to get rid of these before doing analyses of lineage
#' diversification.
#' 
#' Unlike \code{timePaleoPhy}, cal3 methods will always resolve polytomies. In
#' general, this is done using the rate calibrated algorithm, although if
#' argument \code{randres = TRUE}, polytomies will be randomly resolved with uniform
#' probability, ala multi2di from ape. Also, cal3 will always add the terminal
#' ranges of taxa. However, because of the ability to infer potential
#' ancestor-descendant relationships, the length of terminal branches may be
#' shorter than taxon ranges themselves, as budding may have occurred during
#' the range of a morphologically static taxon. By resolving polytomies with
#' the cal3 method, this function allows for taxa to be ancestral to more than
#' one descendant taxon. Thus, users who believe their dataset may contain
#' indirect ancestors are encouraged by the package author to try cal3 methods
#' with their consensus trees, as opposed to using the set of most parsimonious
#' trees. Comparing the results of these two approaches may be very revealing.
#' 
#' Like \code{timePaleoPhy}, \code{cal3TimePaleoPhy} is designed for direct application to datasets 
#' where taxon first and last appearances are precisely known in continuous time, with 
#' no stratigraphic uncertainty. This is an uncommon form of data to have from the fossil record, 
#' although not an impossible form (micropaleontologists often have very precise 
#' range charts, for example). This means that most users \emph{should not} use \code{cal3TimePaleoPhy} directly, 
#' unless they have written their own code to deal with stratigraphic uncertainty. For
#' some groups, the more typical 'first' and 'last' dates represent the minimum
#' and maximum absolute ages for the fossil collections that a taxon is known
#' is known from. Presumably, the first and last appearances of that taxon in
#' the fossil record is at unknown dates within these bounds. These should not
#' be mistaken as the FADs and LADs desired by \code{cal3TimePaleoPhy}, as \code{cal3TimePaleoPhy} 
#' will use the earliest dates provided to calibrate node ages, which is either
#' an overly conservative approach to time-scaling or fairly nonsensical.
#'
#' Alternatively to using \code{cal3TimePaleoPhy}, \code{bin_cal3TimePaleoPhy} is a wrapper of 
#' \code{cal3TimePaleoPhy} which produces time-scaled trees for datasets which only have 
#' interval data available. For each output tree, taxon first and last appearance 
#' dates are placed within their listed intervals under a uniform distribution. 
#' Thus, a large sample of time-scaled trees will approximate the uncertainty in 
#' the actual timing of the FADs and LADs. 
#'
#' The input \code{timeList} object can have overlapping (i.e. non-sequential) intervals,
#' and intervals of uneven size. Taxa alive in the modern should be listed as last 
#' occurring in a time interval that begins at time 0 and ends at time 0. If taxa 
#' occur only in single collections (i.e. their first and last appearance in the 
#' fossil record is synchronous, the argument point.occur will force all taxa
#' to have instantaneous durations in the fossil record. Otherwise, by default,
#' taxa are assumed to first and last appear in the fossil record at different points
#' in time, with some positive duration. The sites matrix can be used to force
#' only a portion of taxa to have simultaneous first and last appearances.
#'
#' By setting the argument \code{nonstoch.bin} to \code{TRUE} in \code{bin_cal3TimePaleoPhy}, the
#' dates are NOT stochastically pulled from uniform bins but instead FADs are
#' assigned to the earliest time of whichever interval they were placed in and
#' LADs are placed at the most recent time in their placed interval. This
#' option may be useful for plotting. The sites argument becomes arbitrary if
#' nonstoch.bin is \code{TRUE}.
#' 
#' If \code{timeData} or the elements of \code{timeList} are actually data frames (as output
#' by \code{read.csv} or \code{read.table}), these will be coerced to a matrix.
#' 
#' A tutorial for applying the time-scaling functions in paleotree,
#' particularly the cal3 method, along with an example using real (graptolite)
#' data, can be found at the following link:
#'
#' http://nemagraptus.blogspot.com/2013/06/a-tutorial-to-cal3-time-scaling-using.html
#' 
#' @rdname cal3TimePaleoPhy
#' @aliases cal3TimePaleoPhy bin_cal3TimePaleoPhy cal3

#' @inheritParams timePaleoPhy

#' @param brRate Either a single estimate of the instanteous rate of branching
#' (also known as the 'per-capita' origination rate, or speciation rate if
#' taxonomic level of interest is species) or a vector of per-taxon estimates
#' @param extRate Either a single estimate of the instanteous extinction rate
#' (also known as the 'per-capita' extinction rate) or a vector of per-taxon
#' estimates
#' @param sampRate Either a single estimate of the instanteous sampling rate or
#' a vector of per-taxon estimates

#' @param ntrees Number of time-scaled trees to output

#' @param anc.wt Weighting against inferring ancestor-descendant relationships.
#' The argument anc.wt allows users to change the default consideration of
#' anc-desc relationships. This value is used as a multiplier applied to the
#' probability of choosing any node position which would infer an
#' ancestor-descendant relationship. By default, anc.wt=1, and thus these
#' probabilities are unaltered. if anc.wt is less than 1, the probabilities
#' decrease and at anc.wt=0, no ancestor-descendant relationships are inferred
#' at all. Can be a single value or a vector of per-taxon values, such as if a
#' user wants to only allow plesiomorphic taxa to be ancestors.

#' @param dateTreatment This argument controls the interpretation of timeData. The default setting
#' 'firstLast' treats the dates in timeData as a column of precise first and last appearances.
#' A second option, added by great demand, is 'minMax' which
#' treats these dates as minimum and maximum bounds on single point dates. Under this option,
#' all taxa in the analysis will be treated as being point dates, such that the first appearance
#' is also the last. These dates will be pulled under a uniform distribution. If 'minMax' is used,
#' add.term becomes meaningless, and the use of it will return an error message. A third option
#' is 'randObs'. This assumes that the dates in the matrix are first and last appearance times,
#' but that the desired time of observation is unknown. Thus, this is much like 'firstLast' except
#' the effective time of observation (the taxon's LAD under 'firstLast') is treated an uncertain date, and is randomly
#' sampled between the first and last appearance times. The FAD still is treated as a fixed number, used
#' for dating the nodes. In previous versions of paleotree, this
#' was called in \code{cal3timePaleoPhy} using the argument rand.obs, which has been removed
#' for clarity. This temporal uncertainty in times of observation might be useful if
#' a user is interested in applying phylogeny-based approaches to studying trait evolution, but have
#' per-taxon measurements of traits that come from museum specimens with uncertain temporal placement.
#' With both arguments 'minMax' and 'randObs', the sampling of dates from random distributions should
#' compel users to produce many time-scaled trees for any given analytical purpose.
#' Note that 'minMax' returns an error in 'bin' time-scaling functions; please use
#' 'points.occur' instead.

# @param rand.obs Should the tips represent observation times uniform
# distributed within taxon ranges? This only impacts the location of tip-dates, 
# i.e. the 'times of observation' for taxa, and does not impact the dates used to 
# determine node ages. Thus, this is an alternative to using only the LADs or only the FADs
# as the per-taxon times of observation. If rand.obs is TRUE, then it is assumed
# that users wish the tips to represent observations made with some temporal
# uncertainty, such that they might have come from any point within a taxon's
# known range.  This might be the case, for example, if a user is interested in
# applying phylogeny-based approaches to studying trait evolution, but have
# per-taxon measurements of traits that come from museum specimens with
# uncertain temporal placement. When rand.obs is TRUE, the tips are placed randomly
# within taxon ranges, as if uniformly distributed, and thus multiple trees should be 
# created and analysed.

#' @param node.mins The minimum dates of internal nodes (clades) on a phylogeny can be set
#' using node.mins. This argument takes a vector of the same length as the number of nodes,
#' with dates given in the same order as nodes are ordered in the \code{tree$edge} matrix.
#' Note that in \code{tree$edge}, terminal tips are given the first set of numbers
#' (\code{1:Ntip(tree)}), so the first element of \code{node.mins} is the first internal node
#' (the node numbered \code{Ntip(tree)+1}, which is generally the root for most \code{phylo}
#' objects read by \code{read.tree}). Not all nodes need be given minimum dates; those without
#' minimum dates can be given as NA in \code{node.mins}, but the vector must be the same length
#' as the number of internal nodes in \code{tree}. These are minimum date constraints, such that
#' a node will be 'frozen' by the cal3 algorithm so that constrained nodes will always be
#' \emph{at least as old as this date}, but the final date may be even older depending on the
#' taxon dates used, the parameters input for the cal3 algorithm and any other minimum node dates
#' given (e.g. if a clade is given a very old minimum date, this will (of course) over-ride any
#' minimum dates given for clades that that node is nested within). if the constrained nodes include
#' a polytomy, this polytomy will still be resolved with respect to the cal3 algorithm, but the first
#' divergence will be 'frozen' so that it is at least as old as the minimum age, while any additional
#' divergences will be allowed to occur after this minimum age.

#' @param FAD.only Should the tips represent observation times at the start of
#' the taxon ranges? If TRUE, result is similar to when terminal ranges are no
#' added on with timePaleoPhy. If FAD.only is TRUE and dateTreatment is 'minMax' or 'randObs', the
#' function will stop and a warning will be produced, as these combinations imply
#' contradictory sets of times of observation.

#' @param adj.obs.wt If the time of observation are before the LAD of a taxon,
#' should the weight of the time of observation be adjusted to account for the
#' known observed history of the taxon which occurs AFTER the time of
#' observation? Should only have an effect if time of observation ISN'T the
#' LAD, if times of observation for a potential ancestor are earlier than the
#' first appearance of their potential descendants and if the ancestral weights
#' for taxa aren't set to zero (so there can be potential ancestors).

#' @param root.max Maximum time before the first FAD that the root can be
#' pushed back to.

#' @param step.size Step size of increments used in zipper algorithm to assign
#' node ages.

#' @param randres Should polytomies be randomly resolved using multi2di in ape
#' rather than using the cal3 algorithm to weight the resolution of polytomies 
#' relative to sampling in the fossil record?

#' @param plot If true, plots the input, "basic" time-scaled phylogeny (an
#' intermediary step in the algorithm) and the output
#' cal3 time-scaled phylogeny.

#' @param tolerance Acceptable amount of shift in tip dates from dates listed
#' in \code{timeData}. Shifts outside of the range of \code{tolerance} will
#' cause a warning message to be issued that terminal tips appear to be
#' improperly aligned.

#' @param diagnosticMode If \code{TRUE}, \code{cal3timePaleoPhy} will return to the
#' console and to graphic devices an enormous number of messages, plots and
#' anciliary information that may be useful or entirely useless to figuring
#' out what is going wrong.

#' @param nonstoch.bin If true, dates are not stochastically pulled from
#' uniform distributions. See below for more details.

#' @return The output of these functions is a time-scaled tree or set of
#' time-scaled trees, of either class phylo or multiphylo, depending on the
#' argument ntrees. All trees are output with an element \code{$root.time}. This is
#' the time of the root on the tree and is important for comparing patterns
#' across trees.
#'
#' Additional elements are \code{sampledLogLike} and \code{$sumLogLike} which respectively
#' record a vector containing
#' the 'log-densities' of the various node-ages selected for each tree by the 'zipper'
#' algorithm, and the sum of those log-densities. Although they are very similar to
#' log-likelihood values, they may not be true likelihoods, as node ages are conditional on the other
#' ages selected by other nodes. However, these values may give an indication about the relative
#' optimality of a set of trees output by the cal3 functions.
#'
#' Trees created with \code{bin_cal3TimePaleoPhy} will output with some additional
#' elements, in particular \code{$ranges.used}, a matrix which records the
#' continuous-time ranges generated for time-scaling each tree. (Essentially a
#' pseudo-timeData matrix.)

#' @note Most importantly, please note the stochastic element of the three
#' rate-calibrated time-scaling methods. These do not use traditional
#' optimization methods, but instead draw divergence times from a distribution
#' defined by the probability of intervals of unobserved evolutionary history.
#' This means analyses MUST be done over many cal3 time-scaled trees for
#' analytical rigor! No one tree is correct.
#'
#' Similarly, please account for stratigraphic uncertainty in your analysis.
#' Unless you have exceptionally resolved data, use a wrapper with the cal3
#' function, either the provided bin_cal3TimePaleoPhy or code a wrapper
#' function of your own that accounts for stratigraphic uncertainty in 
#' your dataset. Remember that the FADs (earliest dates) given to timePaleoPhy
#'  will *always* be used to calibrate node ages!

#' @author David W. Bapst

#' @seealso \code{\link{timePaleoPhy}}, \code{\link{binTimeData}},
#' \code{\link{make_durationFreqCont}}, \code{sProb2sRate}, \code{\link{multi2di}}

#' @references 
#' Bapst, D. W. 2013. A stochastic rate-calibrated method for time-scaling
#' phylogenies of fossil taxa. \emph{Methods in Ecology and Evolution}.
#' 4(8):724-733.
#' 
#' Foote, M. 2000. Origination and extinction components of taxonomic
#' diversity: general problems. Pp. 74-102. In D. H. Erwin, and S. L. Wing,
#' eds. Deep Time: Paleobiology's Perspective. The Paleontological Society,
#' Lawrence, Kansas.
#' 
#' Foote, M. 2001. Inferring temporal patterns of preservation, origination,
#' and extinction from taxonomic survivorship analysis. \emph{Paleobiology}
#' 27(4):602-630.
#' 
#' Friedman, M., and M. D. Brazeau. 2011. Sequences, stratigraphy and
#' scenarios: what can we say about the fossil record of the earliest
#' tetrapods? \emph{Proceedings of the Royal Society B: Biological Sciences}
#' 278(1704):432-439.
#' 
#' Stanley, S. M. 1979. Macroevolution: Patterns and Process. W. H. Freeman,
#' Co., San Francisco.

#' @examples
#' 
#' #Simulate some fossil ranges with simFossilRecord
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #simulate a fossil record with imperfect sampling with sampleRanges
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
#' cladogram <- taxa2cladogram(taxa,plot=TRUE)
#' #this library allows one to use rate calibrated type time-scaling methods (Bapst, in prep.)
#' #to use these, we need an estimate of the sampling rate (we set it to 0.5 above)
#' likFun<-make_durationFreqCont(rangesCont)
#' srRes<-optim(parInit(likFun),likFun,lower=parLower(likFun),upper=parUpper(likFun),
#'       method="L-BFGS-B",control=list(maxit=1000000))
#' sRate <- srRes[[1]][2]
#' # we also need extinction rate and branching rate
#'    # we can get extRate from getSampRateCont too
#' #we'll assume extRate=brRate (ala Foote et al., 1999); may not always be a good assumption
#' divRate<-srRes[[1]][1]
#' #now let's try cal3TimePaleoPhy, which time-scales using a sampling rate to calibrate
#' #This can also resolve polytomies based on sampling rates, with some stochastic decisions
#' ttree <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
#'     sampRate=sRate,ntrees=1,plot=TRUE)
#' #notice the warning it gives!
#' phyloDiv(ttree)
#' 
#' #by default, cal3TimePaleoPhy may predict indirect ancestor-descendant relationships
#' #can turn this off by setting anc.wt=0
#' ttree <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
#'     sampRate=sRate,ntrees=1,anc.wt=0,plot=TRUE)
#' 
#' \donttest{
#' #let's look at how three trees generated with very different time of obs. look
#' ttreeFAD <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
#'     FAD.only=TRUE,dateTreatment="firstLast",sampRate=sRate,ntrees=1,plot=TRUE)
#' ttreeRand <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
#'     FAD.only=FALSE,dateTreatment="randObs",sampRate=sRate,ntrees=1,plot=TRUE)
#' #by default the time of observations are the LADs
#' ttreeLAD <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
#'     FAD.only=FALSE,dateTreatment="randObs",sampRate=sRate,ntrees=1,plot=TRUE)
#' layout(1:3)
#' parOrig <- par(no.readonly=TRUE)
#' par(mar=c(0,0,0,0))
#' plot(ladderize(ttreeFAD));text(5,5,"time.obs=FAD",cex=1.5,pos=4)
#' plot(ladderize(ttreeRand));text(5,5,"time.obs=Random",cex=1.5,pos=4)
#' plot(ladderize(ttreeLAD));text(5,5,"time.obs=LAD",cex=1.5,pos=4)
#' layout(1); par(parOrig)
#' 
#' #to get a fair sample of trees, let's increase ntrees
#' ttrees <- cal3TimePaleoPhy(cladogram,rangesCont,brRate=divRate,extRate=divRate,
#'     sampRate=sRate,ntrees=9,plot=FALSE)
#' #let's compare nine of them at once in a plot
#' layout(matrix(1:9,3,3))
#' parOrig <- par(no.readonly=TRUE)
#' par(mar=c(0,0,0,0))
#' for(i in 1:9){plot(ladderize(ttrees[[i]]),show.tip.label=FALSE)}
#' layout(1)
#' par(parOrig)
#' #they are all a bit different!
#' 
#' #can plot the median diversity curve with multiDiv
#' multiDiv(ttrees)
#' 
#' #using node.mins
#' #let's say we have (molecular??) evidence that node #5 is at least 1200 time-units ago
#' #to use node.mins, first need to drop any unshared taxa
#' droppers <- cladogram$tip.label[is.na(
#'       match(cladogram$tip.label,names(which(!is.na(rangesCont[,1])))))]
#' cladoDrop <- drop.tip(cladogram, droppers)
#' # now make vector same length as number of nodes
#' nodeDates <- rep(NA, Nnode(cladoDrop))
#' nodeDates[5]<-1200
#' ttree <- cal3TimePaleoPhy(cladoDrop,rangesCont,brRate=divRate,extRate=divRate,
#'     sampRate=sRate,ntrees=1,node.mins=nodeDates,plot=TRUE)
#' 
#' #example with time in discrete intervals
#' set.seed(444)
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #simulate a fossil record with imperfect sampling with sampleRanges()
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
#' cladogram <- taxa2cladogram(taxa,plot=TRUE)
#' #Now let's use binTimeData to bin in intervals of 1 time unit
#' rangesDisc <- binTimeData(rangesCont,int.length=1)
#' #we can do something very similar for the discrete time data (can be a bit slow)
#' likFun<-make_durationFreqDisc(rangesDisc)
#' spRes<-optim(parInit(likFun),likFun,lower=parLower(likFun),upper=parUpper(likFun),
#'       method="L-BFGS-B",control=list(maxit=1000000))
#' sProb <- spRes[[1]][2]
#' #but that's the sampling PROBABILITY per bin, not the instantaneous rate of change
#' #we can use sProb2sRate() to get the rate. We'll need to also tell it the int.length
#' sRate1 <- sProb2sRate(sProb,int.length=1)
#' #we also need extinction rate and branching rate (see above)
#'     #need to divide by int.length...
#' divRate<-spRes[[1]][1]/1
#' #estimates that r=0.3... kind of low (simulated sampling rate is 0.5)
#' #Note: for real data, you may need to use an average int.length (no constant length)
#' ttree <- bin_cal3TimePaleoPhy(cladogram,rangesDisc,brRate=divRate,extRate=divRate,
#'     sampRate=sRate1,ntrees=1,plot=TRUE)
#' phyloDiv(ttree)
#' #can also force the appearance timings not to be chosen stochastically
#' ttree1 <- bin_cal3TimePaleoPhy(cladogram,rangesDisc,brRate=divRate,extRate=divRate,
#'     sampRate=sRate1,ntrees=1,nonstoch.bin=TRUE,plot=TRUE)
#' phyloDiv(ttree1)
#'
#' # testing node.mins in bin_cal3TimePaleoPhy
#' ttree <- bin_cal3TimePaleoPhy(cladoDrop,rangesDisc,brRate=divRate,extRate=divRate,
#'     sampRate=sRate1,ntrees=1,node.mins=nodeDates,plot=TRUE)
#' # with randres=TRUE
#' ttree <- bin_cal3TimePaleoPhy(cladoDrop,rangesDisc,brRate=divRate,extRate=divRate,
#'     sampRate=sRate1,ntrees=1,randres=TRUE,node.mins=nodeDates,plot=TRUE)
#' 
#' 
#' #example with multiple values of anc.wt
#' ancWt <- sample(0:1,nrow(rangesDisc[[2]]),replace=TRUE)
#' names(ancWt)<-rownames(rangesDisc[[2]])
#' ttree1 <- bin_cal3TimePaleoPhy(cladogram,rangesDisc,brRate=divRate,extRate=divRate,
#'     sampRate=sRate1,ntrees=1,anc.wt=ancWt,plot=TRUE)
#' }
#' 

#NOT NEEDED
#NULL
# @rdname cal3Timescaling

#' @export
cal3TimePaleoPhy<-function(tree, timeData, brRate, extRate, sampRate,
	ntrees=1, anc.wt=1, node.mins=NULL, dateTreatment="firstLast",
	FAD.only=FALSE, adj.obs.wt=TRUE, root.max=200, step.size=0.1,
	randres=FALSE, noisyDrop=TRUE, tolerance=0.0001, diagnosticMode=FALSE, plot=FALSE){
	
	#function for Ps - use pqr2Ps
	#example data
	#tree<-rtree(10);tree$edge.length<-sample(0:1,Nedge(tree),replace=TRUE);tree<-di2multi(tree)
	#sampRate=rep(0.1,Ntip(tree));names(sampRate)<-tree$tip.label
	#brRate<-extRate<-sampRate
	#timeData<-runif(Ntip(tree),200,400);timeData<-cbind(timeData,timeData-runif(Ntip(tree),1,80))
	#rownames(timeData)<-tree$tip.label
	#ntrees=1; anc.wt=1; node.mins=NULL; dateTreatment="firstLast";
	#FAD.only=FALSE; adj.obs.wt=TRUE; root.max=200; step.size=0.1;
	#randres=FALSE; noisyDrop=TRUE; plot=FALSE
	#
	# record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
	#	nTotalTaxa=c(50,100), nExtant=0)
	# taxa<-fossilRecord2fossilTaxa(record)
	#cladogram<-taxa2cladogram(taxa);timeData<-sampleRanges(taxa,r=0.1)
	#
	###trying to see if adj.wts works
	#tree<-rtree(2);anc.wt=1;node.mins=NULL;brRate<-extRate<-sampRate<-0.1;timeData<-cbind(c(100,95),c(90,85));adj.obs.wt=TRUE
	#rownames(timeData)<-tree$tip.label;root.max=200;plot=TRUE;rand.obs=FALSE;FAD.only=TRUE;ntrees=1;randres=FALSE;step.size=0.1
	#
	#add.zombie=FALSE;node.mins<-c(-sort(-runif(1,600,900)),rep(NA,Nnode(tree)-1))	#assume two very deep divergences
	#
	#require(ape)#;require(phangorn)

	if(!inherits(tree, "phylo")){stop("tree is not of class phylo")}
	if(!inherits(timeData,"matrix")){
		if(inherits(timeData,"data.frame")){
			timeData<-as.matrix(timeData)
		}else{
			stop("timeData not of matrix or data.frame format")
			}
		}
	if(!any(dateTreatment==c("firstLast","minMax","randObs"))){
		stop("dateTreatment must be one of 'firstLast', 'minMax' or 'randObs'!")}
	if(dateTreatment=="randObs" & FAD.only){stop("FAD.only=TRUE and dateTreatment='randObs' are conflicting arguments")}
	if(dateTreatment=="minMax" & FAD.only){
		stop("FAD.only=TRUE and dateTreatment='minMax' are conflicting, as there are no FADs, as dates are simply point occurrences")}
	#first clean out all taxa which are NA or missing in timeData
	if(ntrees==1){message("Warning: Do not interpret a single cal3 time-scaled tree, regardless of other arguments!")}
	if(ntrees<1){stop("ntrees<1")}
	#originalInputTree<-tree
	droppers<-tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeData[,1])))))]
	if(length(droppers)>0){
		if(length(droppers)==Ntip(tree)){
			stop("Absolutely NO valid taxa shared between the tree and temporal data!")}
		if(noisyDrop){
			message(paste("Warning: Following taxa dropped from tree:",
				paste0(droppers,collapse=", ")))}
		tree<-drop.tip(tree,droppers)
		if(Ntip(tree)<2){
			stop("Less than two valid taxa shared between the tree and temporal data!")}
		timeData[which(!sapply(rownames(timeData),function(x) any(x==tree$tip.label))),1]<-NA
		}
	if(!is.null(node.mins)){
		if(length(droppers)>0){	#then... the tree has changed unpredictably, node.mins unusable
			stop("node.mins not compatible with datasets where some taxa are dropped; drop before analysis instead")}
		if(Nnode(tree)!=length(node.mins)){
			stop("node.mins must be same length as number of nodes in the input tree!")}
		}
	timeData<-timeData[!is.na(timeData[,1]),]
	if(any(is.na(timeData))){stop("Weird NAs in Data??")}
	if(any(timeData[,1]<timeData[,2])){stop("timeData is not in time relative to modern (decreasing to present)")}
	if(length(unique(rownames(timeData)))!=length(rownames(timeData))){stop("Duplicate taxa in timeList[[2]]")}
	if(length(unique(tree$tip.label))!=length(tree$tip.label)){stop("Duplicate tip taxon names in tree$tip.label")}
	if(length(rownames(timeData))!=length(tree$tip.label)){
		stop("Odd irreconcilable mismatch between timeList[[2]] and tree$tip.labels")}
	#make sure all taxa have a sampRate, brRate and exRate
	if(length(sampRate)==1){sampRate<-rep(sampRate,Ntip(tree));names(sampRate)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(sampRate))))){
			stop("Sampling Rates Not Given For All Taxa on Tree!")}}
	if(length(brRate)==1){brRate<-rep(brRate,Ntip(tree));names(brRate)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(brRate))))){
			stop("Branching Rates Not Given For All Taxa on Tree!")}}
	if(length(extRate)==1){extRate<-rep(extRate,Ntip(tree));names(extRate)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(extRate))))){
			stop("Extinction Rates Not Given For All Taxa on Tree!")}}
	#allow per-taxon anc.wt values
	if(length(anc.wt)==1){anc.wt<-rep(anc.wt,Ntip(tree));names(anc.wt)<-tree$tip.label
		#if it is a species-named vector, all the species better be there1
		}else{if(any(is.na(match(tree$tip.label,names(anc.wt))))){
			stop("Ancestral Weights Not Given For All Taxa on Tree!")}}
	Ps<-sapply(tree$tip.label,function(x) pqr2Ps(brRate[x],extRate[x],sampRate[x]))
	names(Ps)<-tree$tip.label
	#timescale with timePaleoPhy to get "basic" timetree
	ttree1<-timePaleoPhy(tree=tree,timeData=timeData,
		type="basic",dateTreatment="firstLast",
		node.mins=node.mins,add.term=FALSE,inc.term.adj=FALSE)
	#identify which nodes are min-locked; make sure to update when resolving polytomies
	if(length(node.mins)>0){
		locked_nodesOrig<-which(!is.na(node.mins))+Ntip(tree)
	}else{
		locked_nodesOrig<-NA
		}
	ttree1<-collapse.singles(ttree1)
	ttrees<-rmtree(ntrees,3)
	sampledLogLike<-numeric()
	for(ntr in 1:ntrees){
		#10/30/12: get FAD, new LAD (time of observation), and then calculate difference between t.obs and LAD
		if(dateTreatment=="randObs" | FAD.only){
			if(FAD.only){
				timeData1<-cbind(timeData[,1],timeData[,1],timeData[,1]-timeData[,2])
				}
			if(dateTreatment=="randObs"){
				timeData1<-cbind(timeData[,1],apply(timeData,1,function(x) runif(1,x[2],x[1])))
				timeData1<-cbind(timeData1,timeData1[,2]-timeData[,2])
				}
		}else{
			if(dateTreatment=="minMax"){
				datesUniffy<-apply(timeData,1,function(x) runif(1,x[2],x[1]))
				timeData1<-cbind(datesUniffy,datesUniffy,0)
				}
			timeData1<-cbind(timeData,0)
			}
		# now randomly resolve if randres
		if(randres & (!is.binary.tree(ttree1) | !is.rooted(ttree1))){
			ktree<-multi2di(ttree1)
			# need to updated node locks if any node.mins given
			if(!identical(locked_nodesOrig,NA)){
				origDesc<-lapply(prop.part(ttree1),function(x) sort(ttree1$tip.label[x]))
				treeDesc<-lapply(prop.part(ktree),function(x) sort(ktree$tip.label[x]))
				node_changes<-match(origDesc,treeDesc)
				locked_nodes<-node_changes[locked_nodesOrig]
			}else{
				locked_nodes<-locked_nodesOrig
				}
		}else{
			ktree<-ttree1
			locked_nodes<-locked_nodesOrig
			}
		nodes<-(1:Nnode(ktree))+Ntip(ktree)		#get a vector of all internal nodes	
		nodes<-nodes[order(-node.depth(ktree)[-(1:Ntip(ktree))])]	#order by depth
		anags<-character();budds<-character();nAdjZip<-0
		while(length(nodes)>0){		#can't use a for() because # of nodes may change
			if(diagnosticMode){
				save_tree<-ktree
				dev.new()
				plot(ktree)
				}
			node<-nodes[1]
			tipl<-ktree$tip.label
			#put together tip data: diffLAD is the difference between the time of observation and the true LAD
			tipd<-cbind(ID=(1:Ntip(ttree1)),FAD=(timeData1[tipl,1]),time.obs=(timeData1[tipl,2]),SR=sampRate[tipl],
				BR=brRate[tipl],ER=extRate[tipl],Ps=Ps[tipl],ancWt=anc.wt[tipl],diffLAD=(timeData1[tipl,3]))
			if(node==(Ntip(ktree)+1)){
				min_zip<-(-root.max)	#if root, allow to be push back up to root.max
				stem_len<-root.max
				#root_push<--seq(min_zip,0,by=step.size)
			}else{									#if not root, push down to lower node
				min_zip<-(-ktree$edge.length[ktree$edge[,2]==node])
				stem_len<-ktree$edge.length[ktree$edge[,2]==node]
				}		
			dnodes<-ktree$edge[ktree$edge[,1]==node,2]	#find the daughter nodes
			dlen<-ktree$edge.length[match(dnodes,ktree$edge[,2])]	#find the dedges lengths
			minlocked<-ifelse(!all(is.na(locked_nodes)),any(node==locked_nodes),FALSE)	#is this node min-locked?
			#if(any(dnodes==which(ktree$tip.label=="t10"))){break()}
			if(length(dnodes)>2){		#if node is a polytomy, use PARALLEL ZIPPER
				#pick a starting lineage
				dSR<-drng<-dBR<-dER<-dPs<-danc.wt<-ddfLAD<-numeric()
				for(i in dnodes){		#for each desc, get vector of SR for earliest and range if desc is a tip
					dtips<-match(unlist(Descendants(ktree,i)),tipd[,1])
					dearly<-which(tipd[dtips,2]==max(tipd[dtips,2]))[1]
					dSR[length(dSR)+1]<-tipd[dearly,4]
					dBR[length(dBR)+1]<-tipd[dearly,5]
					dER[length(dER)+1]<-tipd[dearly,6]
					dPs[length(dPs)+1]<-tipd[dearly,7]
					danc.wt[length(danc.wt)+1]<-tipd[dearly,8]
					ddfLAD[length(ddfLAD)+1]<-tipd[dearly,9]		#10-30-12: diff between t.obs and LAD
					drng[length(drng)+1]<-ifelse(length(dtips)>1,NA,diff(unlist(tipd[dtips,3:2])))
					}
				#08-01-12: choice of starting lineage doesn't matter (see SRC method)
				dnode1<-dnodes[which(dlen==min(dlen))[1]]		#just pick first appearing
				#make sure to include stem length in calculations!
				#08-03-12: as with new SRC above, the stem length will JUST be the -min_zip: max root.push!
				#if(node==(Ntip(ktree)+1)){
				#	#07-29-12: this treats the stem length as a single gap
				#		#given that this should be a single gap and not gamma distributed
				#	root_density<-dexp(root_push,rate=dSR[dnode1==dnodes]+(dBR[dnode1==dnodes]*dPs[dnode1==dnodes]))
				#	stem_len<-sample(root_push,1,prob=root_density)
				#	}
				#make data structure for placed lineages; anc= row of anc lineage, events in time-from-stem 
				plin<-c(dnode1,(dlen[dnode1==dnodes]+stem_len),drng[dnode1==dnodes],NA,
					0,dlen[dnode1==dnodes]+stem_len,dlen[dnode1==dnodes]+drng[dnode1==dnodes]+stem_len,
					danc.wt[dnode1==dnodes],ddfLAD[dnode1==dnodes])
				plin<-matrix(plin,1,)
				colnames(plin)<-c("node","brl","rng","anc","tSpec","tFO","tLO","ancWt","diffLAD")
				#place additional lineages, in order of max zip, using parallel zippers along placed lineages
				#add_nodes will hold necessary information on dlen, rng, SR for unplaced nodes
				add_nodes<-cbind(dnodes,dlen+stem_len,drng,dSR,dBR,dER,dPs,danc.wt,ddfLAD)[-which(dnodes==dnode1),]
				add_nodes<-add_nodes[order(dlen[-which(dnodes==dnode1)]),]
				colnames(add_nodes)<-c("node","dstem","rng","SR","BR","ER","Ps","ancWt","diffLAD")	#dstem is distance from stem to FAD
				for(i in 1:nrow(add_nodes)){
					#identify each currently placed lineage, produce zipper for each relative to stem point
					zips<-matrix(,1,3)
					for(j in 1:nrow(plin)){
						min_zip<-plin[j,5]	#time of speciation (stem time for each placed lineage)
						max_zip<-ifelse(plin[j,8]>0 & !is.na(plin[j,7]),	#max zip is complicated, dependent on anc.wt
							min(add_nodes[i,2],plin[j,7]),	#min of LAD of the placed lineage or FAD of the lineage to be placed
							min(add_nodes[i,2],plin[j,6]))	#min of FAD of the placed lineage or FAD of the lineage to be placed
						if(minlocked){max_zip<-stem_len}		#if minlocked, node must be not earlier than original node time
						#and on a stemTime=0 scale as used for the parallel zipper, stem_len is the original node time, unlike single zipper below
						poss_zip<-seq(min_zip,max_zip,by=step.size)
						#08-01-12: getting the density via the Cal3 algorithm
						gap1<-plin[j,6]-poss_zip		#waiting time from FAD1 to zip
						gap1<-ifelse(gap1>0,gap1,0)		#waiting time from branching point to FAD1
						gap2<-add_nodes[i,2]-poss_zip			#waiting time from branching point to FAD2
						gapStem2zip<-ifelse(plin[j,6]>poss_zip,poss_zip,plin[j,6])	#waiting time from stem to FAD1 OR br node
						totalgap<-gap1+gap2+gapStem2zip
						#07-31-12: Given a lack of other options, gamma(shape=2,rate=r+p*Ps) distribution best fit
							#under different combinations with p=q=r,p=q>r and p=q<r
							#admittedly not a perfect fit, though!
							#use rates from the lineages being added
						linDense<-dgamma(totalgap,shape=2,rate=add_nodes[i,4]+(add_nodes[i,5]*add_nodes[i,7]))
						linDense<-ifelse(poss_zip>plin[j,6],
							plin[j,8]*linDense,	#with anc.wt (of the PLACED LINEAGE)
								    linDense)	#without anc.wt, but with gap1 probability
						#10-30-12 need to correct linDense value for time of observation if difference from LAD =/= 0
							#use the diffLAD of the PLACED LINEAGE (CAUSE ITS THE ANCESTOR)
							#if adj.obs.wt, if anc.wt>0 & diffLAD>0 & plin's LAD is before the add_nodes's FAD...
						if(ifelse(is.na(plin[j,7]),FALSE,add_nodes[i,2]>(plin[j,7]+step.size)) & 
								adj.obs.wt & plin[j,8]>0 & plin[j,9]>0){	
							adj_max<-max(add_nodes[i,2],plin[j,9]+plin[j,7])
							adj_zips<-seq(plin[j,7]+step.size,adj_max,by=step.size) 	#the zips we ain't looking at
							#08-01-12: getting the density via the Cal3 algorithm
							gap1<-plin[j,6]-adj_zips		#waiting time from FAD1 to zip
							gap1<-ifelse(gap1>0,gap1,0)		#waiting time from branching point to FAD1
							gap2<-add_nodes[i,2]-adj_zips			#waiting time from branching point to FAD2
							gapStem2zip<-ifelse(plin[j,6]>adj_zips,adj_zips,plin[j,6])	#waiting time from stem to FAD1 OR br node
							adj_totalgap<-gap1+gap2+gapStem2zip
							#07-31-12: Given a lack of other options, gamma(shape=2,rate=r+p*Ps) distribution best fit
								#under different combinations with p=q=r,p=q>r and p=q<r
								#admittedly not a perfect fit, though!
							adj_linDense<-dgamma(adj_totalgap,shape=2,rate=add_nodes[i,4]+(add_nodes[i,5]*add_nodes[i,7]))
							linDense[length(linDense)]<-linDense[length(linDense)]+sum(adj_linDense*plin[j,8])
							nAdjZip<-nAdjZip+1
							}
						new_zip<-cbind(linDense,plin[j,1],poss_zip)
						zips<-rbind(zips,new_zip)
						}
					if(nrow(zips)<3){zips<-matrix(zips[-1,],1,3)}else{zips<-zips[-1,]}
					colnames(zips)<-c("linDensity","anc","tzip")
					#new as of 08-21-12
					zip_prob<-zips[,1]
					zip_prob[is.na(zip_prob)]<-0
					if(sum(zip_prob)==0){zip_prob<-rep(1,length(zip_prob))}
					zip_prob<-zip_prob/sum(zip_prob)
					ch_zip<-sample(1:nrow(zips),1,prob=zip_prob)	#sample zips
					#record sampled zip_prob
					sampledLogLike<-c(sampledLogLike,log(zip_prob[ch_zip]))
					ch_anc<-zips[ch_zip,2]
					ch_tzip<-zips[ch_zip,3]
					#if anagenesis, add to anags; if budding, add to budds
					if(!is.na(plin[ch_anc==plin[,1],7])){	#if the anc is terminal
						if(plin[ch_anc==plin[,1],7]==ch_tzip){anags<-c(anags,ktree$tip.label[ch_anc])}	#if anagenetic
						if(plin[ch_anc==plin[,1],6]<ch_tzip){budds<-c(budds,ktree$tip.label[ch_anc])}		#if budding
						}
					new_lin<-c(add_nodes[i,1],add_nodes[i,2]-ch_tzip,add_nodes[i,3],
						ch_anc,ch_tzip,add_nodes[i,2],add_nodes[i,2]+add_nodes[i,3],add_nodes[i,8],add_nodes[i,9])
					plin<-rbind(plin,new_lin)	#put in plin
					}
				#turn into a subtree using taxa2phylo()
				taxad_o<-t(apply(plin,1,function(x) c(x[1],x[4],x[5],ifelse(is.na(x[7]),x[6],x[7]))))
				new_anc<-sapply(taxad_o[,2],function(x) ifelse(is.na(x),NA,which(taxad_o[,1]==x)))
				taxad_n<-cbind(1:nrow(taxad_o),new_anc,taxad_o[,3:4])
				taxad_n[,3:4]<-max(taxad_n[,3:4])-taxad_n[,3:4]
				rownames(taxad_n)<-paste("t",1:nrow(taxad_o),sep="")
				subtree<-taxa2phylo(taxad_n)
				subtree$tip.label<-taxad_o[match(subtree$tip.label,paste("t",1:nrow(taxad_o),sep="")),1]
				new_stem<-diff(sort(plin[,5]))[1]	#time to stem (branch length for stem)
				#stick desc tips onto the subtree
				for(i in dnodes){
					dtip<-which(subtree$tip.label==i)
					if(i>Ntip(ktree)){		#if its a clade
						subclade<-extract.clade(ktree,i)
						subtree<-bind.tree(subtree,subclade,where=dtip)
						subtree<-collapse.singles(subtree)
					}else{				#if its a tip
						subtree$tip.label[dtip]<-ktree$tip.label[i]
					}}
				#replace original node with new, resolved, scaled node
				if(node!=(Ntip(ktree)+1)){	#if it isn't the node
					drtips<-prop.part(ktree)[[node-Ntip(ktree)]]
					tip_lab<-ktree$tip.label[drtips[1]]	#I need to cut out all but one tip, for the sole purpose of putting it all back later)
					droptree<-collapse.singles(drop.tip(ktree,drtips[-1]))
					droptree$edge.length[droptree$edge[,2]==which(droptree$tip.label==tip_lab)]<-new_stem	#reset edge length leading to remaining tip to new_stem
					droptree<-bind.tree(droptree,subtree,where=which(droptree$tip.label==tip_lab))	#put in subtree at tip
					ktree1<-droptree
				}else{				#if it is the node
					ktree1<-subtree
					}				
				#once you've changed the structure of the tree find original nodes in new tree (surprisingly frustating to code!)
				if(length(nodes)>1){
					d_o<-lapply(Descendants(ktree,nodes[-1]),function(x) ktree$tip.label[x])
					d_n<-lapply(Descendants(ktree1)[-(1:Ntip(ktree1))],function(x) ktree1$tip.label[x])
					nodes1<-sapply(d_o,function(x) which(sapply(d_n,function(y) 
						ifelse(length(y)==length(x),all(sort(y)==sort(x)),FALSE))))
					nodes1<-Ntip(ktree1)+nodes1
					nodes1<-nodes1[order(-node.depth(ktree1)[nodes1])]	#order by depth
					if(!all(is.na(locked_nodes))){	#update locked_nodes, can re-use d_n
						d_ol<-lapply(Descendants(ktree,locked_nodes),function(x) ktree$tip.label[x])
						locked_nodes<-sapply(d_ol,function(x) which(sapply(d_n,function(y) 
							ifelse(length(y)==length(x),all(sort(y)==sort(x)),FALSE))))
						locked_nodes<-Ntip(ktree1)+locked_nodes
						}					
				}else{nodes1<-numeric()}				#don't bother if no more nodes left...
				if(diagnosticMode){
					layout(matrix(1:2,2,));plot(save_tree);plot(ktree1);layout(1)
					}
				#update tipd and nodes (tree str will have changed)
				ktree1<-collapse.singles(ktree1)
				ktree<-ktree1
				nodes<-nodes1
			}else{	#if node is NOT a polytomy, then use regular zipper
				dlen1<-min(dlen)				#the shortest branch (nothing can be done about this one)
				dlen2<-max(dlen)				#the longest branch
				dnode1<-dnodes[which(dlen==dlen1)[1]]
				dnode2<-dnodes[dnodes!=dnode1]
				#get the SR of earliest tip of d2
				d2FADs<-tipd[match(unlist(Descendants(ktree,dnode2)),tipd[,1]),2]	#need the NODE for prop.part, idiot!
				d2early<-which(d2FADs==max(d2FADs))	
				d2early<-ifelse(length(d2early)>1,sample(d2early,1),d2early)	#if more than one of same FAD, just randomly choose one
				d2SR<-(tipd[match(unlist(Descendants(ktree,dnode2)),tipd[,1]),4])[d2early]
				d2BR<-(tipd[match(unlist(Descendants(ktree,dnode2)),tipd[,1]),5])[d2early]
				d2Ps<-(tipd[match(unlist(Descendants(ktree,dnode2)),tipd[,1]),7])[d2early]
				d1rng<-ifelse(dnode1<=Ntip(ktree),diff(unlist(tipd[match(dnode1,tipd[,1]),3:2])),NA)	#if clade, range = NA	
				d2rng<-ifelse(dnode2<=Ntip(ktree),diff(unlist(tipd[match(dnode2,tipd[,1]),3:2])),NA)
				d1ancWt<-ifelse(dnode1<=Ntip(ktree),unlist(tipd[match(dnode1,tipd[,1]),8]),0)
				d1diffLAD<-ifelse(dnode1<=Ntip(ktree),unlist(tipd[match(dnode1,tipd[,1]),9]),0)
				#ZIPPPER: first create a list of scenarios, treat the position of the node like a zipper which can be moved up or down
				max_zip<-ifelse(d1ancWt>0 & !is.na(d1rng),min(dlen1+d1rng,dlen2),dlen1)	#if d1 isn't a clade and there can be ancestors	
				if(minlocked){max_zip<-0}		#in single zipper, unlike parallel zipper, 0 is original node time
				poss_zip<-seq(min_zip,max_zip,by=step.size)
				#waiting time from zip (branching point) to FAD1
				gap1<-ifelse(poss_zip>dlen1,0,dlen1-poss_zip)		#restrict to be dlen1 if longer than dlen1 (FAD1-time, prob 0 anyway!)
				gap2<-dlen2-poss_zip			#waiting time from branching point to FAD2
				#waiting time from stem to FAD1 OR branch node (zip)
				gapStem2zip<-ifelse((stem_len+dlen1)>(stem_len+poss_zip),
					stem_len+poss_zip,stem_len+dlen1)	
				totalgap<-gap1+gap2+gapStem2zip
				#07-31-12: Given a lack of other options, gamma(shape=2,rate=r+p*Ps) distribution best fit
					#under different combinations with p=q=r,p=q>r and p=q<r
						#admittedly not a perfect fit, though!
				linDense<-dgamma(totalgap,shape=2,rate=d2SR+(d2BR*d2Ps))
				linDensity<-ifelse((stem_len+dlen1)<(stem_len+poss_zip),
					d1ancWt*linDense,	#with anc.wt (of first taxon)
					 	  linDense)	#without anc.wt, but with gap1 probability
				#10-30-12 need to correct linDense value for time of observation if difference from LAD =/= 0
					#use the diffLAD of the potential ancestor
					#if adj.obs.wt, if anc.wt>0 & diffLAD>0 & d1's LAD is before d2's FAD...
				if(adj.obs.wt & d1ancWt>0 & d1diffLAD>0 & (dlen1+d1rng+step.size)<dlen2){
					adj_max<-max(dlen2,d1diffLAD+(dlen1+d1rng))
					adj_zips<-seq((dlen1+d1rng)+step.size,adj_max,by=step.size) 	#the zips we ain't looking at
					#08-01-12: getting the density via the Cal3 algorithm
					#waiting time from zip (branching point) to FAD1
					gap1<-ifelse(adj_zips>dlen1,0,dlen1-adj_zips)		#restrict to be dlen1 if longer than dlen1 (FAD1-time, prob 0 anyway!)
					gap2<-dlen2-adj_zips			#waiting time from branching point to FAD2
					#waiting time from stem to FAD1 OR branch node (zip)
					gapStem2zip<-ifelse((stem_len+dlen1)>(stem_len+adj_zips),
						stem_len+adj_zips,stem_len+dlen1)
					adj_totalgap<-gap1+gap2+gapStem2zip
					#07-31-12: Given a lack of other options, gamma(shape=2,rate=r+p*Ps) distribution best fit
						#under different combinations with p=q=r,p=q>r and p=q<r
						#admittedly not a perfect fit, though!
					adj_linDense<-dgamma(totalgap,shape=2,rate=d2SR+(d2BR*d2Ps))
					linDensity[length(linDensity)]<-linDensity[length(linDensity)]+sum(adj_linDense*d1ancWt)
					nAdjZip<-nAdjZip+1
					}
				#new as of 08-21-12
				linDensity[is.na(linDensity)]<-0
				if(sum(linDensity)==0){linDensity[length(linDensity)]<-1}
				linDensity1<-linDensity/sum(linDensity)
				ch_zip<-sample(poss_zip,1,prob=linDensity1)	#pick zipper location
				#record sampled zip_prob
				sampledLogLike<-c(sampledLogLike,log(linDensity1[ch_zip]))
				#calculate new branch lengths, adding terminal ranges to tips
				new_dlen1<-ifelse(ch_zip>dlen1,NA,dlen1-ch_zip)			#If not budding or anagenesis
				new_dlen2<-dlen2-ch_zip		
				if(is.na(new_dlen1)){							#if budding or anagensis
					if(ch_zip==max_zip & dlen1+d1rng<dlen2){ 			#if anagensis
						anags<-c(anags,ktree$tip.label[dnode1])
						new_dlen1<-0
					}else{								#if budding
						budds<-c(budds,ktree$tip.label[dnode1])
						new_dlen1<-d1rng+dlen1-ch_zip
				}}else{new_dlen1<-ifelse(!is.na(d1rng),new_dlen1+d1rng,new_dlen1)}	#if tip, add term range
				new_dlen2<-ifelse(!is.na(d2rng),new_dlen2+d2rng,new_dlen2)		#if tip, add rng to new dlen
				#rescale branches according to their new lengths
				if(node!=(Ntip(ktree)+1)){	#if not root
					ktree$edge.length[ktree$edge[,2]==node]<-ch_zip-min_zip	#change stem length
					}
				ktree$edge.length[match(dnode1,ktree$edge[,2])]<-new_dlen1
				ktree$edge.length[match(dnode2,ktree$edge[,2])]<-new_dlen2
				if(diagnosticMode){
					print(c(node,ch_zip-min_zip,new_dlen1,new_dlen2))
					}
				nodes<-nodes[-1]			#update nodes
			}}
		ktree<-reorder(collapse.singles(ktree),"cladewise")
		ktree$anag.tips<-anags	#record the number of anagenetic ancestors
		ktree$budd.tips<-budds	#record the number of budding ancestors	
		ktree$nAdjZip<-nAdjZip
		#record sampled log-likelihoods
		ktree$sampledLogLike<-sampledLogLike
		ktree$sumLogLike<-sum(sampledLogLike)
		#now add root.time: because NO TIPS ARE DROPPED (due to anagenesis) can calculate this now
			#must be calculated on LADs because the terminal ranges are added to the TREE!!!
			#should be time of earliest LAD + distance of root from earliest tip
		ktree$root.time<-max(timeData1[ktree$tip.label,2])+min(node.depth.edgelength(ktree)[1:Ntip(ktree)])	
		names(ktree$edge.length)<-names(ktree$tip.label)<-names(ktree$budd.tips)<-names(ktree$anag.tips)<-NULL
		#stuff for checking if things are correct
		tipdiffs<-cbind(	diff(sort(-timeData1[,2])),
					diff(sort(node.depth.edgelength(ktree)[1:Ntip(ktree)])),
					diff(sort(-timeData1[,2]))-diff(sort(node.depth.edgelength(ktree)[1:Ntip(ktree)])))	
		test1<-all(tipdiffs[,3]<tolerance)
		test2<-identical(names(sort(-timeData1[,2])),
				ktree$tip.label[order(node.depth.edgelength(ktree)[1:Ntip(ktree)])])
		if(length(unique(timeData1[,2]))<Ntip(tree)){
			test2<-TRUE}	#test 2 does not work if any LADS are same
		if(all(c(test1,test2))){
			ktree$test<-"passed"
		}else{
			warning("Warning: Terminal tips improperly aligned, cause unknown. Use ouput with care.")
			}
		if(plot){
			parOrig <- par(no.readonly=TRUE)
			par(mar=c(2.5,2.5,1,2.5))
			layout(matrix(1:3,3,))
			plot(ladderize(tree),show.tip.label=TRUE,use.edge.length=FALSE)
			plot(ladderize(ttree1),show.tip.label=TRUE);axisPhylo()			
			plot(ladderize(ktree),show.tip.label=TRUE);axisPhylo()
			layout(1)
			par(parOrig)		
			}
		names(ktree$edge.length)<-NULL
		ttrees[[ntr]]<-ktree
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}

#' @rdname cal3TimePaleoPhy
#' @export
bin_cal3TimePaleoPhy<-function(tree,timeList,brRate,extRate,sampRate,
	ntrees=1, anc.wt=1, node.mins=NULL, dateTreatment="firstLast",
	FAD.only=FALSE,sites=NULL,point.occur=FALSE,nonstoch.bin=FALSE,
	adj.obs.wt=TRUE,root.max=200,step.size=0.1,
	randres=FALSE,noisyDrop=TRUE,
	tolerance=0.0001, diagnosticMode=FALSE, plot=FALSE){
	#see the bin_cal3 function for more notation...
	#require(ape)
	if(!inherits(tree, "phylo")){
		stop("tree is not of class phylo")
		}
	if(!inherits(timeList[[1]],"matrix")){
		if(inherits(timeList[[1]],"data.frame")){
			timeList[[1]]<-as.matrix(timeList[[1]])
		}else{
			stop("timeList[[1]] not of matrix or data.frame format")
			}
			}
	if(!inherits(timeList[[2]],"matrix")){
		if(inherits(timeList[[2]],"data.frame")){
			timeList[[2]]<-as.matrix(timeList[[2]])
		}else{
			stop("timeList[[2]] not of matrix or data.frame format")
			}
		}
	if(dateTreatment=="minMax"){stop("Instead of dateTreatment='minMax', please use argument point.occur instead in bin_ functions or use sites argument")}
	if(!any(dateTreatment==c("firstLast","randObs"))){
		stop("dateTreatment must be one of 'firstLast' or 'randObs'!")}
	if(ntrees<1){stop("ntrees<1")}
	if(ntrees==1){message("Warning: Do not interpret a single cal3 time-scaled tree")}
	if(ntrees==1 & !nonstoch.bin){
		message("Warning: Do not interpret a single tree; dates are stochastically pulled from uniform distributions")}
	if(dateTreatment=="randObs" & FAD.only){stop("FAD.only=TRUE and dateTreatment='randObs' are conflicting arguments")}
	if(dateTreatment=="minMax" & FAD.only){
		stop("FAD.only=TRUE and dateTreatment='minMax' are conflicting, as there are no FADs, as dates are simply point occurrences")}
	if(!is.null(sites) & point.occur){stop("Inconsistent arguments, point.occur=TRUE would replace input 'sites' matrix\n Why not just make site assignments for first and last appearance the same in your input site matrix?")}
	#clean out all taxa which are NA or missing for timeData
	droppers<-tree$tip.label[is.na(match(tree$tip.label,names(which(!is.na(timeList[[2]][,1])))))]
	if(length(droppers)>0){
		if(length(droppers)==Ntip(tree)){stop("Absolutely NO valid taxa shared between the tree and temporal data!")}
		if(noisyDrop){
			message(paste("Warning: Following taxa dropped from tree:",paste0(droppers,collapse=", ")))
			}
		tree<-drop.tip(tree,droppers)
		if(is.null(tree)){stop("Absolutely NO valid taxa shared between the tree and temporal data!")}
		if(Ntip(tree)<2){stop("Less than two valid taxa shared between the tree and temporal data!")}
		timeList[[2]][which(!sapply(rownames(timeList[[2]]),function(x) any(x==tree$tip.label))),1]<-NA
		}
	if(!is.null(node.mins)){
		if(length(droppers)>0){	#then... the tree has changed unpredictably, node.mins unusable
			stop("node.mins not compatible with datasets where some taxa are dropped; drop before analysis instead")}
		if(Nnode(tree)!=length(node.mins)){
			stop("node.mins must be same length as number of nodes in the input tree!")}
		}
	#best to drop taxa from timeList that aren't represented on the tree
	notTree<-rownames(timeList[[2]])[is.na(match(rownames(timeList[[2]]),tree$tip.label))]
	if(length(notTree)>0){
		if(is.null(sites)){
			if(noisyDrop){
				message(paste("Warning: Following taxa dropped from timeList:",paste0(notTree,collapse=", ")))}
			timeList[[2]]<-timeList[[2]][!is.na(match(rownames(timeList[[2]]),tree$tip.label)),]
		}else{
			stop("Some taxa in timeList not included on tree: no automatic taxon drop if 'sites' are given. Please remove from both sites and timeList and try again.")
			}
		}	
	timeList[[2]]<-timeList[[2]][!is.na(timeList[[2]][,1]),]
	if(any(is.na(timeList[[2]]))){
		stop("Weird NAs in Data??")}
	if(any(apply(timeList[[1]],1,diff)>0)){
		stop("timeList[[1]] not in intervals in time relative to modern")}
	if(any(timeList[[1]][,2]<0)){
		stop("Some dates in timeList[[1]] <0 ?")}
	if(any(apply(timeList[[2]],1,diff)<0)){
		stop("timeList[[2]] not in intervals numbered from first to last (1 to infinity)")}
	if(any(timeList[[2]][,2]<0)){
		stop("Some dates in timeList[[2]] <0 ?")}
	if(length(unique(rownames(timeList[[2]])))!=length(rownames(timeList[[2]]))){
		stop("Duplicate taxa in timeList[[2]]")}
	if(length(unique(tree$tip.label))!=length(tree$tip.label)){
		stop("Duplicate tip taxon names in tree$tip.label")}
	if(length(rownames(timeList[[2]]))!=length(tree$tip.label)){
		stop("Odd irreconcilable mismatch between timeList[[2]] and tree$tip.labels")}
	if(is.null(sites)){
		if(point.occur){
			if(any(timeList[[2]][,1]!=timeList[[2]][,2])){
				stop("point.occur=TRUE but some taxa have FADs and LADs listed in different intervals?!")}
			sites<-matrix(c(1:Ntip(tree),1:Ntip(tree)),Ntip(tree),2)
		}else{
			sites<-matrix(1:(nrow(timeList[[2]])*2),,2)
			}
	}else{	#make sites a bunch of nicely behaved sorted integers
		sites[,1]<-sapply(sites[,1],function(x) which(x==sort(unique(as.vector(sites)))))
		sites[,2]<-sapply(sites[,2],function(x) which(x==sort(unique(as.vector(sites)))))
		}
	ttrees<-rmtree(ntrees,3)
	siteTime<-matrix(,max(sites),2)
	for (i in unique(as.vector(sites))){		#build two-col matrix of site's FADs and LADs
		go<-timeList[[2]][which(sites==i)[1]]	#find an interval for this site
		siteTime[i,]<-timeList[[1]][go,]
		}
	for(ntrb in 1:ntrees){
		if(!nonstoch.bin){
			bad_sites<-unique(as.vector(sites))
			siteDates<-apply(siteTime,1,function(x) runif(1,x[2],x[1]))
			while(length(bad_sites)>0){
				siteDates[bad_sites]<-apply(siteTime[bad_sites,],1,function(x) runif(1,x[2],x[1]))
				bad_sites<-unique(as.vector(sites[(siteDates[sites[,1]]-siteDates[sites[,2]])<0,]))
				#message(length(bad_sites))
				}
			timeData<-cbind(siteDates[sites[,1]],siteDates[sites[,2]])
		}else{
			timeData<-cbind(siteTime[sites[,1],1],siteTime[sites[,2],2])
			}
		rownames(timeData)<-rownames(timeList[[2]])
		#if(rand.obs){timeData[,2]<-apply(timeData,1,function(x) runif(1,x[2],x[1]))}
		#if(FAD.only){timeData[,2]<-timeData[,1]}
		tree2<-suppressMessages(
			cal3TimePaleoPhy(tree,timeData,brRate=brRate,extRate=extRate,sampRate=sampRate,
				ntrees=1,anc.wt=anc.wt,node.mins=node.mins,adj.obs.wt=adj.obs.wt,
				root.max=root.max,step.size=step.size,
				FAD.only=FAD.only,dateTreatment=dateTreatment,randres=randres,
				tolerance=tolerance,diagnosticMode=diagnosticMode,plot=plot)
			)
		colnames(timeData)<-c("FAD","LAD")
		tree2$ranges.used<-timeData
		names(tree2$edge.length)<-NULL
		ttrees[[ntrb]]<-tree2
		}
	if(ntrees==1){ttrees<-ttrees[[1]]}
	return(ttrees)
	}