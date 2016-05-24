#' Simulation of sequencing read counts and genotypes
#'
#' Simulates genotypes and read counts under the model of Blischak \emph{et al}.
#'
#' Total reads are simulated using a Poisson distribution with mean equal to the \code{coverage} set by the user.
#' Next, genotypes are simulated for the specified number of individuals using the vector of allele frequencies provided to the function.
#' The number of loci simulated is equal to the number of elements supplied by the vector of allele frequencies.
#' The number of reference reads is then simulated using Eq. 1 from Blischak \emph{et al}. using the total reads, genotypes and sequencing error.
#'
#' @param pVec A vector of allele frequencies strung together using the concatenate function.
#' @param N_ind The number of individuals to simulate.
#' @param coverage The average number of sequences simulated per individual per locus (Poisson distributed).
#' @param ploidy The ploidy level of individuals in the population.
#' @param error The level of sequencing error. A fixed constant.
#'
#' @return A list of 3 matrices:
#' \describe{
#'  \item{\code{genos}}{A matrix of the simulated genotypes.}
#'  \item{\code{tot_read_mat}}{A matrix of the simulated number of total reads.}
#'  \item{\code{ref_read_mat}}{A matrix of the simulated number of reference reads.}
#'  }
#'
#' @references Blischak PD, Kubatko LS, Wolfe AD. 2015. Accounting for genotype uncertainty in the estimation of allele frequencies in autopolyploids. \emph{In review}. bioRxiv, \strong{doi}:####.
#' @useDynLib polyfreqs
#' @importFrom Rcpp sourceCpp

#' @export
sim_reads <- function(pVec, N_ind, coverage, ploidy, error){
  genos <- matrix(apply(as.matrix(pVec), 1, function(x) rbinom(N_ind, ploidy, x)), nrow=N_ind, ncol=length(pVec))
  tot_read_mat <- matrix(rpois(N_ind*length(pVec), coverage),nrow=N_ind, ncol=length(pVec))
  ref_read_mat <- sim_ref_reads(tot_read_mat, genos, ploidy, error)

  rownames(tot_read_mat) <- paste("ind", 1:N_ind, sep="")
  rownames(ref_read_mat) <- paste("ind", 1:N_ind, sep="")

  colnames(tot_read_mat) <- paste("loc", 1:length(pVec), sep="")
  colnames(ref_read_mat) <- paste("loc", 1:length(pVec), sep="")

  return(list(genos=genos, tot_read_mat=tot_read_mat, ref_read_mat=ref_read_mat))
}
