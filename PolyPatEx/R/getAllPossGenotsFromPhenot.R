
##  Function getAllPossGenotsFromPhenot
##
## Given a phenotype and a ploidy, generate all possible genotypes.
##
## \code{getAllPossGenotsFromPhenot} generates all possible genotypes
## from a specified phenotype and ploidy.
##
## This routine uses \code{\link{combinations}} from the \pkg{gtools}
## package.
##
## @title Generate the genotypes from a phenotype
## @param phenot character vector a phenotype, as a vector of
## alleles.
## @param p integer the ploidy of the species, an even number greater
## than or equal to 4.
## @return A list of character vectors, each of length equal to the
## ploidy, and each specifying a possible genotype arising from the
## original phenotype.
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## \dontrun{
##
## getAllPossGenotsFromPhenot(phenot=c("P1","M1"),p=4)
##
## }
##
getAllPossGenotsFromPhenot <- function(phenot,p) {
  ##
  n <- length(phenot)
  if (n==p) {  ## The phenotype IS a genotype
    return(list(phenot))
  } else {
    ## If the phenotype is not also a genotype, then there are
    ##  duplicated alleles in the genotype.  Obtain all possible
    ##  combo's of these duplicated alleles.
    ## [
    ##  Aside: combinations(n,r,v,repeats) : Obtain all distinct
    ##          (unordered) draws of r elements, with/without repeats,
    ##          from a vector v containing n (unique) objects
    ##          (n = length(v))
    ## ]
    dupAlleles <- gtools::combinations(n=n, r=p-n, v=phenot,
                                       repeats.allowed=TRUE)
    dupAList <- split(dupAlleles, row(dupAlleles)) ##Convert to a list
    names(dupAList) <- NULL
    ## Combine the duplicate possibilities with the alleles in the
    ## phenotype, return the sorted allelesets
    allGenotypes <- lapply(dupAList,
                           function (vv,phenot) {
                             sort(c(phenot,vv))},
                           phenot)
    return(allGenotypes)
  }
}
