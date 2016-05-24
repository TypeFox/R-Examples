#' @title Mean difference of a discrete random variable
#' @description Compute the Gini mean difference of a discrete random variable.
#' @usage deltaprob(m, prob)
#' @aliases deltaprob
#' @param m Number of categories
#' @param prob Probability distribution of the random variable
#' @keywords univar
#' @export deltaprob
#' @examples
#' m<-7
#' prob<-c(0.04,0.04,0.05,0.10,0.21,0.32,0.24)
#' deltaprob(m,prob)


deltaprob <-
function(m,prob){
  frip<-cumsum(prob)
  frip1<-frip[1:(m-1)]
  2*sum(frip1*(1-frip1))
}
