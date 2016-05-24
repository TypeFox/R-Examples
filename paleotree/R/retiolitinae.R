#' Cladogram and Range Data for the Retiolitinae
#' 
#' The majority rule consensus cladogram for 22 genera from the Retiolitinae, a
#' clade of Silurian retiolitids, along with discrete time interval data
#' taken from the same publication (Bates et al., 2005). Additional character state
#' data are included for three major, binary-state morphological traits.

#' 
#' @details Interval dates were taken from Sadler et al. (2009). These zones were not a
#' 1-1 match to those in Bates et al., so it took some merging and splitting by
#' the package author, so buyer beware.
#'
#' Character data are from an in prep manuscript containing character data for certain
#' major morphological innovations of graptoloids, coded for a large number of genera based
#' on an extensive survey of the published descriptions. The character data presented here
#' is a small subset of the full dataset.

#' @name retiolitinae

#' @rdname retiolitinae

#' @aliases retiolitinae retioRanges retioTree retioChar

#' @docType data

#' @format This dataset is composed of three objects:
#'
#' \describe{

#'  \item{\code{retioTree}}{An \code{ape} 'phylo' object containing the consensus cladogram.} 

#'  \item{\code{retioRanges}}{A list containing two matrices. The first matrix describes the first
#' and last interval times for 20 Silurian graptolite zones and the second matrix describes when the
#' various genera on the cladogram first and last appear in those graptolite zones. (In other words,
#' \code{retioRanges} has the 'timeList' format called by some paleotree functions).}

#' \item{\code{retioChar}}{A matrix containing binary presence-absence character states for these 22
#' Retiolitinae genera for three characters which they vary in: the presence of extrathecal threads
#' (note only one taxon lacks this character state), the presence of determinant growth and the
#' secondary loss of a nema via resorbtion. Note these character do not vary within these genera.}

#'}

#' @source 
#' Source for cladogram and zonal ranges for genera: 
#'
#' Bates, D. E. B., A. Kozlowska, and A. C. Lenz. 2005. Silurian retiolitid graptolites:
#' Morphology and evolution. \emph{Acta Palaeontologica Polonica} 50(4):705-720.
#' 
#' Source for interval dates for graptolite zones: 
#'
#' Sadler, P. M., R. A. Cooper, and M. Melchin. 2009. High-resolution, early Paleozoic (Ordovician-Silurian)
#' time scales. \emph{Geological Society of America Bulletin} 121(5-6):887-906.
#'
#' Source for morphological character data:
#'
#' Collected for Bapst and Mitchell, in prep. 

#' @seealso For more example graptolite datasets, see \code{\link{graptDisparity}}

#' @keywords datasets

#' @examples
#' 
#' #load data
#' data(retiolitinae)
#' 
#' #Can plot discrete time interval diversity curve with retioRanges
#' taxicDivDisc(retioRanges)
#' 
#' #Can plot the unscaled cladogram
#' plot(retioTree)
#' #Can plot the determinant growth character on the cladogram
#' tiplabels(pch=16,col=(retioChar[,2]+1),adj=0.25)
#' 
#' #Use basic time-scaling (terminal branches only go to FADs)
#' ttree<-bin_timePaleoPhy(tree=retioTree,timeList=retioRanges,type="basic",
#' 	ntrees=1,plot=TRUE)
#' 
#' #Note that this function creates stochastic time-scaled trees...
#' 	#A sample of 1 is not representative!
#' 
#' #phylogenetic diversity curve
#' phyloDiv(ttree)
#' 
#' 
NULL