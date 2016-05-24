
##  Function getPossMPPhenotCounts
##
## Generate all possible mother and progeny phenotype counts from a
## specified ploidy.
##
## Given a ploidy \eqn{p = 4, 6, 8, \ldots}, this function generates
## all possible combinations of phenotype counts for mother and
## progeny.
##
## The counts generated are the number of alleles in the mother
## (\code{nM}), the number of alleles in the progeny (\code{nP}) and
## the number of alleles that appear in both mother and progeny
## (\code{nMP}).
##
## The table excludes the trivial cases where \code{nM}, \code{nP} or
## \code{nMP} are zero.  Also note that \code{nMP} will be no greater
## than the minimum of \code{nM} and \code{nP}.
##
## @title Generate possible mother and progeny phenotype counts
## @param p integer: the ploidy, an even number greater than or equal
## to 4.
## @return A data frame containing the columns:
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
## }
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## \dontrun{
##
## getPossMPPhenotCounts(p=4)
##
## }
##
getPossMPPhenotCounts <- function(p) {
  ##
  if ((p < 4) || (p %% 2 > 0.5)) {
    ## Error check : p should be even, and 4 or greater
    stop("\n Ploidy should be an even number, 4 or greater \n")
  }
  ## Possible numbers of alleles in Mother, Progeny
  ddtmp <- expand.grid(nM=1:p,nP=1:p)
  ## Maximum possible shared alleles
  ss <- apply(ddtmp,1,min)
  return(with(ddtmp,data.frame(nM = rep(nM,times=ss),
                               nP = rep(nP,times=ss),
                               nMP = sequence(ss))))
}
