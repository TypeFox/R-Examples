#' @title Locus Summaries
#' @description Compile standard by-locus summaries.
#' 
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata logical. If \code{TRUE}, return a list of summary matrices 
#'   for each stratum.
#' @param ... arguments to be passed on to summary functions.
#' 
#' @return A matrix with rows for each locus and columns containing summaries of:
#'   \tabular{ll}{
#'     \code{num.genotyped} \tab The number of samples genotyped.\cr
#'     \code{prop.genotyped} \tab The proportion of samples genotyped.\cr
#'     \code{num.alleles} \tab The number of alleles in the locus.\cr
#'     \code{allelic.richness} \tab The allelic richness of the locus.\cr
#'     \code{prop.unique.alleles} \tab Proportion of alleles found in a single sample.\cr
#'     \code{expt.heterozygosity} \tab Expected heterozygosity.\cr
#'     \code{obsvd.heterozygosity} \tab Observed heterozygosity.\cr
#'   }
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' data(msats.g)
#' msats.g <- stratify(msats.g, "fine")
#' 
#' summarizeLoci(msats.g)
#' 
#' @export
#' 
summarizeLoci <- function(g, by.strata = FALSE, ...) {
  summary.stats <- function(x) {
    n.genotyped <- nInd(x) - numMissing(x)
    cbind(
      num.genotyped = n.genotyped,
      prop.genotyped = n.genotyped / nInd(x),
      num.alleles = numAlleles(x),
      allelic.richness = allelicRichness(x),
      prop.unique.alleles = propUniqueAlleles(x),
      exptd.heterozygosity = exptdHet(x),
      obsvd.heterozygosity = obsvdHet(x)
    )
  }
  
  if(by.strata) {
    lapply(strataSplit(g), summary.stats)
  } else summary.stats(g)
}