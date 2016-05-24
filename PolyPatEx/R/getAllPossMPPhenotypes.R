
##  Function getAllPossMPPhenotypes
##
## Generate all of the nontrivial mother and progeny phenotypes
## corresponding to the given ploidy.
##
## Given a ploidy (4, 6, 8, \ldots), generates all possible
## non-trivial mother and progeny phenotypes.
##
## \sQuote{Non-trivial} means that cases where there are no maternal
## alleles, or no progeny alleles, or no alleles shared between
## mother and progeny, are ignored.
##
## This function calls \code{\link{getPossMPPhenotCounts}} to get the
## table of possible (non-trivial) phenotype counts, then calls
## \code{\link{getMPPhenotypeFromCount}} on each set of counts to
## generate the corresponding maternal and progeny phenotypes.
## These are returned in a list (the default) or in a tabular form as
## a data frame.
##
## In the data frame form, maternal and progeny allele vectors are
## pasted together to form single character strings, with the alleles
## separated by spaces.
##
## @title Generate all mother and progeny phenotypes
## @param p integer: the ploidy, an even number greater than or equal to 4.
## @param asList logical: return the result as a list (the default), or as a
## data frame (\code{asList=FALSE})?
## @return Either a list or a data frame.
##
## The list form of the output contains an element for each possible
## combination of the phenotype counts, nM, nP & nMP.   Each element
## is itself a list, containing the elements
## \describe{
## \item{nM}{the number of alleles in the maternal phenotype.}
## \item{nP}{the number of alleles in the progeny's phenotype.}
## \item{nMP}{the number of alleles that appear in both maternal and
##  progeny phenotypes.}
## \item{MPhenot}{The maternal phenotype, as a character vector of
##                alleles.}
## \item{PPhenot}{The progeny phenotype, as a character vector of
##                alleles.}
## }
##
## The data frame form of the output has the above components as
## columns in a data frame, where the phenotypes are represented as
## character strings with spaces separating the alleles.
##
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## \dontrun{
##
## getAllPossMPPhenotypes(p=4)
##
## getAllPossMPPhenotypes(p=4,asList=FALSE)
##
## }
##
getAllPossMPPhenotypes <- function(p,asList=TRUE) {
  ##
  ## Get data frame of possible values for nM, nP, nMP
  dd <- getPossMPPhenotCounts(p)
  ## Convert to a list
  paramsList <- apply(dd,1,
                      function(vv) list(nM = vv[1],
                                        nP = vv[2],
                                        nMP = vv[3]))
  ## Obtain all mother and progeny phenotypes
  alleleSetsList <- lapply(paramsList,
                           function(ll) {
                             do.call(getMPPhenotypeFromCount,ll)
                           })
  ## Combine the parameter values with the corresponding phenotypes in
  ##  a single list structure
  AllPossMPPhenotypes <- mapply(
                                function(params,alleleSets) {
                                  c(params,alleleSets)
                                },
                                paramsList,alleleSetsList,SIMPLIFY=FALSE)
  ## If list output is not desired, convert to a data frame
  if (!asList) {
    AllPossMPPhenotypes <- as.data.frame(
                               t(sapply(AllPossMPPhenotypes,
                                        function(ll){
                                          with(ll,c(nM,nP,nMP,
                                                    paste(MPhenot,collapse=" "),
                                                    paste(PPhenot,collapse=" ")))
                                        })),
                                         stringsAsFactors=FALSE)
    names(AllPossMPPhenotypes) <-  c("nM","nP","nMP","MPhenot","PPhenot")
    AllPossMPPhenotypes$nM <- as.numeric(AllPossMPPhenotypes$nM)
    AllPossMPPhenotypes$nP <- as.numeric(AllPossMPPhenotypes$nP)
    AllPossMPPhenotypes$nMP <- as.numeric(AllPossMPPhenotypes$nMP)
  }
  return(AllPossMPPhenotypes)
}
