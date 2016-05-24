
##  Package PolyPatEx
##
##' Polyploid paternity exclusion by allele matching.
##'
##' PolyPatEx provides functions to perform paternity exclusion
##' analysis in autopolyploid species having ploidy 4, 6, or 8.  The
##' package requires codominant marker data from two or more loci in a
##' monoecious or dioecious species.  The marker data can be
##' 'genotypic' (copy numbers known) or 'phenotypic' (copy numbers not
##' known).  PolyPatEx can also perform exclusion on diploid (ploidy
##' 2) \emph{genotype} data.
##'
##' Routines are provided to compare each progeny with its mother, and
##' then with the candidate fathers, to determine which candidates are
##' indeed capable of being fathers, on the basis of the allele sets
##' they display at each locus.
##'
##' PolyPatEx addresses the question - at a given locus, can the
##' candidate father provide a viable gamete given its allele set, and
##' given the possible paternal gametes indicated by the progeny's and
##' mother's allele sets?
##'
##' Note that PolyPatEx does not implement a probabilistic solution to
##' the exclusion problem, merely a simple comparative analysis based
##' on available alleles and their multiplicities.  Also note that
##' PolyPatEx is not optimised for very large marker datasets such as
##' SNP datasets, instead is suited to low density, high information
##' markers such as SSRs.
##'
##' @name PolyPatEx-package
##' @docType package
##' @author Alexander Zwart (alec.zwart at csiro.au)
##' @importFrom gtools combinations
##'
NULL


