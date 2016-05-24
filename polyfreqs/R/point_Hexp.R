#' Estimation of expected heterozygosity
#'
#' \emph{INTERNAL}: Estimates a posterior distribution for the per locus expected heterozygosity using the unbiased estimator of Hardy (2015) and the poterior samples of allele frequencies calculated by \code{\link{polyfreqs}}.
#'
#' Posterior distributions for the per locus expected heterozygosity are automatically calculated and returned by the \code{\link{polyfreqs}} function.
#'
#' @param p_samp A posterior sample of allele frequencies from \code{\link{polyfreqs}}.
#' @param genotypes Matrix of genotypes sampled during MCMC.
#' @param ploidy The ploidy level of individuals in the population (must be >= 2).
#' @return Returns the per locus estimates of expected heterozygosity (\code{per_locus_Hexp})
#' @references Hardy, OJ. 2015. Population genetics of autopolyploids under a mixed mating model and the estimation of selfing rate. \emph{Molecular Ecology Resources}, doi: 10.1111/1755-0998.12431.

#' @export
point_Hexp <- function(p_samp, genotypes, ploidy){

  #Get the number of individuals
  Nind <- nrow(genotypes)


  # Get expected probability of drawing 2 alleles that are the same
  allele1_hom <- p_samp^2
  allele2_hom <- (1 - p_samp)^2

  per_locus_sum_mat <- allele1_hom + allele2_hom

  # Get the observed heterozygosity correction from Hardy (2015)

  obs_heterozygosity <- function(x, ploidy){
    tmp <- (na.omit(x) * (ploidy - na.omit(x)))/choose(ploidy,2)
    return(sum(tmp))
  }

  het_obs <- apply(genotypes, 2, obs_heterozygosity, ploidy=ploidy)
  obs_het_sum <- (ploidy-1)/(ploidy*(Nind^2))*het_obs

  # Calculate the per locus expected heterozygosity with the correction for small sample size
  per_locus_Hexp <- (Nind/(Nind-1))*(1 - per_locus_sum_mat - obs_het_sum)

  return(per_locus_Hexp)
}
