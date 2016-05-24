#' @title Probability distribution of a CUSH model
#' @aliases probcush
#' @description Compute the probability distribution of a CUSH model without covariates.
#' @keywords distribution
#' @export probcush 
#' @usage probcush(m, delta, shelter)
#' @param m Number of ordinal categories
#' @param delta Shelter parameter
#' @param shelter Category corresponding to the shelter choice
#' @return The vector of the probability distribution of a CUSH model without covariates
#' @seealso \code{\link{CUSH}}
#' @references 
#' Capecchi S. and Piccolo D. (2015). Dealing with heterogeneity/uncertainty in sample survey with 
#' ordinal data, \emph{IFCS Proceedings, University of Bologna} \cr
#' Capecchi S and Iannario M. (2015). Gini heterogeneity index for detecting uncertainty in ordinal data surveys,
#'  \emph{SIS Proceedings, Treviso - Ca'Foscari, University of Venice}
#' @examples
#' m<-10
#' shelter<-1
#' delta<-0.4
#' pr<-probcush(m, delta, shelter)
#' plot(1:m,pr,type="h",xlab="Number of categories")
#' points(1:m,pr,pch=19)



probcush <-
function(m,delta,shelter){
  delta*(ifelse(seq(1,m)==shelter,1,0) - 1/m) + 1/m
}
