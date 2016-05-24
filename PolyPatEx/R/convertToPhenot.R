
##  Function convertToPhenot
##
##' Convert a 'genotypic' dataset (marker dosages known) to
##' a 'phenotypic' dataset (marker dosages not known, unique alleles
##' appear only once in each allele set).
##'
##' In the terminology used by the PolyPatEx package, a 'genotypic'
##' allele dataset is one where marker dosages are known, hence any
##' locus at which fewer than p (the ploidy) alleles are detected is
##' incomplete (and subsequently ignored by the genotype-specific
##' routines in this package). A 'phenotypic' allele dataset is one
##' where marker dosages are not known, hence individual alleles
##' appear only once in an allele set, and a complete allele set can
##' contain between \eqn{1} and \eqn{p} alleles.
##'
##' \code{convertToPhenot} converts a genotypic dataset to a
##' phenotypic dataset, simply by removing any allele duplicates from
##' each allele set.  This is probably not something many will want to
##' do, since one loses considerable information in the process...
##'
##' @title Convert a genotype allele dataset to a phenotype dataset
##' @param adata data frame: a genotypic allele dataset.
##' @return A data frame, containing the phenotypic form of the
##' original genotypic dataset.
##' @author Alexander Zwart (alec.zwart at csiro.au)
##' @export
##' @examples
##'
##' ## Using the example dataset 'FR_Genotype':
##' data(FR_Genotype)
##'
##' ## Since we did not load this dataset using inputData(), we must
##' ## first process it with preprocessData() before doing anything
##' ## else:
##' gData <- preprocessData(FR_Genotype,
##'                         numLoci=7,
##'                         ploidy=4,
##'                         dataType="genotype",
##'                         dioecious=TRUE,
##'                         mothersOnly=TRUE)
##'
##' head(gData)  ## Checked and Cleaned version of FR_Genotype
##'
##' pData <- convertToPhenot(gData)
##'
##' head(pData)
##'
convertToPhenot <- function(adata) {
  ##
  checkForValidPPEDataset(adata)
##  dataType <- attr(adata,"dataType")
  if (attr(adata,"dataType") != "genotype") stop("\n Data is not genotypic\n\n")
  numLoci <- attr(adata,"numLoci")
  ploidy <- attr(adata,"ploidy")
  dioecious <- attr(adata,"dioecious")
##  selfCompatible <- attr(adata,"selfCompatible")
  cat("\n Converting...\n")
  for (locus in 1:numLoci) {
    locusRange <- (3+dioecious) + (locus-1)*ploidy + 1:ploidy
    ##Note: apply() returns the _transpose_ of the desired matrix
    ## below - hence the use of t()...
    adata[ , locusRange] <- t(apply(adata[ , locusRange],
                                     1,
                                     function(vv) {
                                       uu <- stripNAs(unique(vv))
                                         is.na(vv) <- TRUE
                                       if(length(uu)>0) {vv[1:length(uu)] <- sort(uu)}
                                       return(vv)
                                     }))
  }
  attr(adata,"dataType") <- "phenotype"
  cat("\n Done...\n")
  return(adata)
}
