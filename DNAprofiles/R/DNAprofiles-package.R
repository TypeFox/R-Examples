#' Methods for studying forensic DNA profiles, in particular familial databases searches.
#' 
#' Main feautures include sampling of DNA profiles/databases as well as relatives of profiles, database comparison exercises, likelihood ratio computations and there are/will be methods to efficiently evaluate the distribution of likelihood ratios. 
#' 
#' @name DNAprofiles-package
#' @aliases DNAprofiles
#' @docType package
#' @title DNA profiles analysis
#' @author \email{m.v.kruijver@@vu.nl}
#' @keywords package
NULL
#' @name freqsNLngm
#' @title Allelic frequencies at 15 NGM loci 
#' @description Allelic frequencies from a reference population (consisting of 2,085 Dutch males) used by the Netherlands Forensic Institute.
#' @docType data
#' @usage freqsNLngm
#' @format A list containing a vector (loci) and a sub-list (freqs).
#' @source Netherlands Forensic Institute
#' @seealso \code{\link{freqsNLsgmplus}},\code{\link{freqsUScaucs}}
#'
#' 
NULL
#' @name freqsNLsgmplus
#' @title Allelic frequencies at 10 SGMplus loci 
#' @description Allelic frequencies from a reference population (consisting of 2,085 Dutch males) used by the Netherlands Forensic Institute.
#' @docType data
#' @usage freqsNLngm
#' @format A list containing a vector (loci) and a sub-list (freqs).
#' @source Netherlands Forensic Institute
#' @seealso \code{\link{freqsNLngm}},\code{\link{freqsUScaucs}}
#'
#' 
NULL
#' @name freqsUScaucs
#' @title Allelic frequencies at 13 CODIS loci 
#' @description Allelic frequencies from a US reference population. This dataset is copied from the \code{relSim} package, raw data is available from \url{http://www.fbi.gov/about-us/lab/forensic-science-communications/fsc/july1999/dnaloci.txt}
#' @docType data
#' @usage freqsNLngm
#' @format A list containing a vector (loci) and a sub-list (freqs).
#' @source Copied from the relSim package
#' @seealso \code{\link{freqsNLsgmplus}},\code{\link{freqsNLngm}}
#'
#' 
NULL
