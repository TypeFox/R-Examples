#' paleotree: Paleontological and Phylogenetic Analyses of Evolution
#' 
#' Analyzes, time-scales and simulates phylogenies of extinct/fossil lineages,
#' along with calculation of diversity curves. Also fits likelihood models to
#' estimate sampling rates from stratigraphic ranges.
#' 
#' \tabular{ll}{ Package: \tab paleotree\cr Type: \tab Package\cr
#' License: \tab CC0\cr } This package
#' contains functions for analyzing sampling rates given ranges of fossil taxa,
#' in both continuous and discrete time, functions for a posteriori time-scaling phylogenies
#' of fossil taxa and functions for simulating the fossil record in both taxic
#' and phylogenetic varieties.
#' 
#' @name paleotree-package
#' @aliases paleotree-package paleotree
#' @docType package
#' @author David W. Bapst
#' 
#' Maintainer: David W. Bapst <dwbapst@@gmail.com>

#' @seealso This package relies extensively on the phylogenetic toolkit and
#' standards offered by the \code{\link[ape:ape-package]{ape}} package, and hence
#' lists this package as a depends, so it is loaded simultaneously.

#' @references 
#' Bapst, D.W. 2012. paleotree: an R package for paleontological
#' and phylogenetic analyses of evolution. \emph{Methods in Ecology and
#' Evolution}. 3: 803-807. doi: 10.1111/j.2041-210X.2012.00223.x
#' 
#' Bapst, D. W. 2013. A stochastic rate-calibrated method for time-scaling
#' phylogenies of fossil taxa. \emph{Methods in Ecology and Evolution}.
#' 4(8):724-733.
#' 
#' Bapst, D. W. 2013. When Can Clades Be Potentially Resolved with Morphology?
#' \emph{PLoS ONE}. 8(4):e62312.
#'
#' Bapst, D. W. 2014. Assessing the effect of time-scaling methods on
#' phylogeny-based analyses in the fossil record. \bold{Paleobiology}
#' \bold{40}(3):331-351.
#'
#' @examples
#' 
#' # get the package version of paleotree
#' packageVersion("paleotree")
#'
#' # get the citation for paleotree
#' citation("paleotree")
#'
#' ##Simulate some fossil ranges with simFossilRecord
#' set.seed(444);
#' record<-simFossilRecord(p=0.1, q=0.1, nruns=1,
#'	nTotalTaxa=c(30,40), nExtant=0)
#' taxa<-fossilRecord2fossilTaxa(record)
#' #let's see what the 'true' diversity curve looks like in this case
#' #plot the FADs and LADs with taxicDivCont()
#' taxicDivCont(taxa)
#' 
#' #simulate a fossil record with imperfect sampling with sampleRanges
#' rangesCont <- sampleRanges(taxa,r=0.5)
#' #plot the diversity curve based on the sampled ranges
#' layout(1:2)
#' taxicDivCont(rangesCont)
#' #Now let's use binTimeData to bin in intervals of 10 time units
#' rangesDisc <- binTimeData(rangesCont,int.length=10)
#' #plot with taxicDivDisc
#' taxicDivDisc(rangesDisc)
#' #compare to the continuous time diversity curve
#' 
#' layout(1)
#'
#' #taxa2phylo assumes we know speciation events perfectly... what if we don't?
#' #first, let's use taxa2cladogram to get the 'ideal' cladogram of the taxa
#' cladogram <- taxa2cladogram(taxa,plot=TRUE)
#' #Now let's try timePaleoPhy using the continuous range data
#' ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",plot=TRUE)
#' #plot diversity curve
#' phyloDiv(ttree,drop.ZLB=TRUE)
#' 
#' #that tree lacked the terminal parts of ranges (tips stops at the taxon FADs)
#' #let's add those terminal ranges back on with add.term
#' ttree <- timePaleoPhy(cladogram,rangesCont,type="basic",add.term=TRUE,plot=TRUE)
#' #plot diversity curve 
#' phyloDiv(ttree)
#' 

# NAMESPACE IMPORTING

#' @import ape
#' @import stats

#' @importFrom phangorn Descendants Ancestors phyDat ancestral.pars
#' @importFrom phytools bind.tip
#' @importFrom graphics par layout plot hist lines legend polygon title axis
#' @importFrom grDevices rainbow dev.new
#' @importFrom methods is
#' @importFrom utils read.csv type.convert


NULL
