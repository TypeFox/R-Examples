#' Point estimation of allele frequencies based on read counts
#'
#' @description \code{simple_freqs} estimates allele frequencies based on read count ratios.
#' @author Paul Blischak
#' @param tM Total reads matrix: matrix containing the total number of reads mapping to each locus for each individual.
#' @param rM Reference reads marix: matrix containing the number of reference reads mapping to each locus for each individual.
#' @return A vector of allele frequencies, one for each locus. Named \code{allele_freqs_hat}.

#' @export
simple_freqs <- function(tM, rM){
  missing_data<-(tM==0)
  tM[missing_data]<-NA
  ratio_mat <- rM/tM
  allele_freqs_hat <- apply(ratio_mat, 2, function(x) mean(na.omit(x)))

  return(allele_freqs_hat)
}
