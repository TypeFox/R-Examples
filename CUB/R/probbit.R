#' @title Probability distribution of a shifted Binomial random variable
#' @description Return the shifted Binomial probability distribution.
#' @aliases probbit
#' @usage probbit(m, csi)
#' @param m Number of ordinal categories
#' @param csi Feeling parameter 
#' @import stats
#' @export probbit
#' @return The vector of the probability distribution of a shifted Binomial model
#' @keywords distribution
#' @examples
#' m<-7
#' csi<-0.7
#' pr<-probbit(m, csi)
#' plot(1:m, pr, type="h", main="Shifted Binomial probability distribution",xlab="Categories")
#' points(1:m,pr,pch=19)

probbit <-
function(m,csi){dbinom(0:(m-1),m-1,1-csi)}
