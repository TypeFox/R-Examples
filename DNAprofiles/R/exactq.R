#' Compute exact LR exceedance probabilities
#' 
#' @param t numeric (vector), threshold
#' @param dists list of per-locus probability distributions of a likelihood ratio
#' @return numeric (vector) with estimated probabilities
#' @details For a combined likelihood ratio \deqn{LR=LR_1 LR_2 \times LR_m,}
#' define \eqn{q_{t|H}} as the probability that the LR exceeds \eqn{t} under hypothesis \eqn{H}, i.e.:
#' \deqn{q_{t|H} := P(LR>t|H).}
#' The hypothesis \eqn{H} can be \eqn{H_p}, \eqn{H_d} or even another hypothesis. The current function computes \eqn{q_{t|H}} by brute force.
#' 
#' @examples
#'
#' data(freqsNLsgmplus)
#'
#'x <- sample.profiles(N = 1, freqsNLsgmplus)
#'
#'# dist of PI for true parent/offspring pairs
#'hp <- ki.dist(x = x, hyp.1="PO",hyp.2="UN",hyp.true="PO",freqs.ki=freqsNLsgmplus)
#'
#'# dist of PI for unrelated pairs
#'hd <- ki.dist(x = x, hyp.1="PO",hyp.2="UN",hyp.true="UN",freqs.ki=freqsNLsgmplus)
#'
#'set.seed(100)
#'
#'# estimate P(PI>1e4) for true PO
#'sim.q(t=1e4,dists=hp)
#'
#'# estimate P(PI>1e4) for unrelated pairs
#'sim.q(t=1e4,dists=hd) # small probability, so no samples exceed t=1e6
#'
#'# importance sampling can estimate the small probability reliably
#'# by sampling from H_p and weighting the samples appropriately
#'sim.q(t=1e4,dists=hd,dists.sample=hp)
#'
#'# compare to exact values
#'exact.q(t = 1e4, dists=hp)
#'exact.q(t = 1e4, dists=hd)
#'
#' @export
exact.q <- function(t,dists){
  dists <- lapply(dists,check.dist) # check if dists are properly specified
  tmp <- Zpdists.properties(dists)
  sapply(t,Zexactq,x=Zdiststomatrix.X(dists=dists),prob=Zdiststomatrix.P(dists=dists),i=tmp$i0,n=tmp$n0,pr0=tmp$prod.pr.0)
}