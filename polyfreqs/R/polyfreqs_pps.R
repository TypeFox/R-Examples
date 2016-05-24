#' Posterior predictive model checks for polyfreqs
#'
#' Uses the posterior distribution of allele frequences from a \code{\link{polyfreqs}} run to test model fit using the posterior predictive model checking procedure described in Blischak \emph{et al}.
#'
#' The observed read count ratio (r/t) for each locus is summed across individuals and then compared to a distribution of read ratios simulated using the posterior allele frequencies by taking their difference.
#' The criterion for passing/failing the posterior predictive check is then made on a per locus basis based on whether or not the distribution of read ratio differences contains 0 in the 95% higherst posterior density interval.
#'
#' @param p_post A matrix containing the posterior samples from a \code{\link{polyfreqs}} run.
#' @param tM Total reads matrix: matrix containing the total number of reads mapping to each locus for each individual.
#' @param rM Reference reads marix: matrix containing the number of reference reads mapping to each locus for each individual.
#' @param ploidy Ploidy level of individuals in the population.
#' @param error The level of sequencing error. A fixed constant.
#'
#' @return A list with two items:
#' \describe{
#'  \item{ratio_diff}{The posterior predictive samples of the difference between the simulated read ratios and the observed read ratio summed across individuals at each locus.}
#'  \item{locus_fit}{A logical vector indicating whether or not each locus passed or failed the posterior predictive model check.}
#' }
#'
#' @references Blischak PD, LS Kubatko and AD Wolfe. Accounting for genotype uncertainty in the estimation of allele frequencies in autopolyploids. \emph{In revision}.
#'
#' @useDynLib polyfreqs
#' @importFrom Rcpp sourceCpp
#'
#' @export
polyfreqs_pps <- function(p_post, tM, rM, ploidy, error){

  sim_ref_read_ratios <- matrix(NA, nrow=nrow(p_post), ncol=ncol(p_post))
  sim_genos <- matrix(NA, nrow=nrow(tM),ncol=ncol(tM))
  sim_ref_read <- matrix(NA, nrow=nrow(tM),ncol=ncol(tM))
  missing_data<-(tM==0)
  obs_ref_read_ratio <- apply(rM/tM, 2, sum, na.rm=T)

  for(i in 1:nrow(p_post)){
    sim_genos <- matrix(apply(as.matrix(p_post[i,]), 1, function(x) rbinom(nrow(tM), ploidy, x)), nrow=nrow(tM), ncol=ncol(tM))
    sim_genos[missing_data]<-NA
    sim_ref_read <- sim_ref_reads(tM, sim_genos, ploidy, error)
    sim_ref_read_ratios[i,] <- apply(sim_ref_read/tM, 2, sum, na.rm=T)
  }

  ratio_diff <- apply(sim_ref_read_ratios, 1, function(x) (obs_ref_read_ratio - x))
  locus_vec <- apply(t(ratio_diff), 2, function(x) sum(sign(quantile(x,c(0.025,0.975))))==0)
  return(list(ratio_diff=t(ratio_diff), locus_fit=locus_vec))
}
