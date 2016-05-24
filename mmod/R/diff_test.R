#' An exact test of population differentiation for genind objects
#'
#' This function uses Fisher's exact test to determine if alleles in 
#' sub-populations are drawn randomly from a larger population (i.e. a 
#' significance test for allelic differentiation among sub-populations).
#'
#' Note, this test returns p-values for each locus in a dataset _not_  
#' estimates of effect size. Since most populations have some degree of 
#' population differentiation, very large samples are almost guaranteed to 
#' return significant results. Refer to estimates of the various differentiation
#' statistics (D, G''ST and Phi'ST)to ascertain  how meaningful such results 
#' might be.
#'
#' @param x a genind object (from package adegenet)
#' @param sim boolean: if TRUE simulate p-value by using an MCMC sample of 
#' those tables that have the same marginal totals as the observed data 
#' (required for all but the smallest datasets)
#' @param nreps number of steps used to simulate p-value (default 2000)
#' @return named vector of p-values testing the null hypothesis these samples
#' where drawn from a panmictic population.
#' @importFrom stats fisher.test
#' @importFrom stats complete.cases
#' @seealso \code{\link{fisher.test}}, which this function wraps
#' @export
#' @examples
#' 
#' data(nancycats)
#' diff_test(seploc(nancycats)[[2]], nreps=100)
#' 


diff_test <- function(x, sim=TRUE, nreps=2000){
  # The test to be applied to each locus
  per.locus <- function(locus){
    allele_counts <- locus@tab*2
    if(any(is.na(locus@tab))){
      warning(paste("Ommiting individuals with NAs for locus", locus$loc.names))
      include <- complete.cases(allele_counts)
      allele_counts <- allele_counts[include,]
      pops <- pop(locus)[include, drop=TRUE]
    }
    else {
      pops <- pop(locus)
    }
    # Sum up alleles (rows in matrix) by population using tapply
    pop_alleles <- apply(allele_counts, 2, function(a) tapply(a, pops, sum))
    return(fisher.test(pop_alleles, simulate.p.value=sim, B=nreps)$p.value)
  }

  loci <- (sapply(seploc(x), per.locus))
  return(loci)
}
