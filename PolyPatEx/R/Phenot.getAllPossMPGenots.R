
##  Function Phenot.getAllPossMPGenots
##
## Given a ploidy, generate all non-trivial mother and progeny
## genotypes.
##
## Given a ploidy \eqn{p = 4, 6, 8, \ldots}, generate all non-trivial
## maternal and progeny genotypes.
##
## Here, \sQuote{non-trivial} means that cases where there are no
## matching alleles between mother and progeny are ignored.
##
## The generated genotypes are returned in a list structure along
## with the phenotype counts and phenotypes that generated the
## genotypes.
##
## @title Generate all mother and progeny genotypes
## @param p the ploidy, an even number greater than or equal to 4.
## @return A list, with as many elements as there are possible
## mother-progeny phenotype combinations for the given ploidy.
##
## Each element of the list is itself a list, containing the following:
##
## \describe{
##
## \item{\code{nM}}{the number of alleles in the maternal phenotype.}
##
## \item{\code{nP}}{the number of alleles in the progeny's
## phenotype.}
##
## \item{\code{nMP}}{the number of alleles that appear in both
## maternal and progeny phenotypes.}
##
## \item{\code{MPhenot}}{The maternal phenotype, as a character
## vector of alleles.}
##
## \item{\code{PPhenot}}{The progeny phenotype, as a character vector
## of alleles.}
##
## \item{\code{MGenots}}{a list of character vectors, each character
## vector representing a maternal genotype that can be generated from
## the maternal phenotype.}
##
## \item{\code{PGenots}}{a list of character vectors, each character
## vector representing a progeny genotype that can be generated from
## the progeny phenotype.}
##
## }
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## \dontrun{
##
## Phenot.getAllPossMPGenots(p=4)
##
## }
##
Phenot.getAllPossMPGenots <- function(p) {
  ## (Phenot prefix to avoid confusion between dataType cases)
  ##
  ll <- getAllPossMPPhenotypes(p)
  for (i in 1:length(ll)) {
    ll[[i]]$MGenots <- getAllPossGenotsFromPhenot(ll[[i]]$MPhenot,p)
    ll[[i]]$PGenots <- getAllPossGenotsFromPhenot(ll[[i]]$PPhenot,p)
  }
  return(ll)
}
