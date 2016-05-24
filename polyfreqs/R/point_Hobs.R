#' Estimation of observed heterozygosity
#'
#' \emph{INTERNAL}: Estimates a posterior distribution for the per locus observed heterozygosity using the unbiased estimator of Hardy (2015) and the poterior samples of genotypes calculated by \code{\link{polyfreqs}}.
#'
#' Posterior distributions for the per locus observed heterozygosity are automatically calculated and returned by the \code{\link{polyfreqs}} function.
#'
#' @param genotypes A matrix of estimated genotypes returned from the function \code{\link{get_map_genotypes}}.
#' @param ploidy The ploidy level of individuals in the population (must be >= 2).
#' @return Returns per locus estimates of observed heterozygosity (\code{per_locus_Hobs}).
#' @references Hardy, OJ. 2015. Population genetics of autopolyploids under a mixed mating model and the estimation of selfing rate. \emph{Molecular Ecology Resources}, doi: 10.1111/1755-0998.12431.

#' @export
point_Hobs <- function(genotypes, ploidy){

  mean_heterozygosity <- function(x, ploidy){
    tmp <- (na.omit(x) * (ploidy - na.omit(x)))/choose(ploidy,2)
    return(mean(tmp))
  }

  per_locus_Hobs <- apply(genotypes, 2, mean_heterozygosity, ploidy=ploidy)

  return(per_locus_Hobs)
}
