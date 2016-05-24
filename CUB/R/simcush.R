#' @title Simulation routine for CUSH models
#' @aliases simcush
#' @description Generate \eqn{n} pseudo-random numbers following the distribution of a CUSH
#'  model without covariates.
#' @keywords distribution
#' @usage simcush(n, m, delta, shelter)
#' @export simcush
#' @import stats
#' @param n Number of simulated observations
#' @param m Number of ordinal categories
#' @param delta Shelter parameter
#' @param shelter Category corresponding to the shelter choice
#' @seealso  \code{\link{probcush}}
#' @examples
#' n<-200
#' m<-7
#' delta<-0.3
#' shelter<-3
#' simulation<-simcush(n, m, delta, shelter)
#' plot(table(simulation), xlab="Ordinal categories", ylab="Frequencies")


simcush <-
function(n,m,delta,shelter){
  dicodelta<-runif(n)<delta
  uncert<-sample(m,n,replace=TRUE)
  return(dicodelta*shelter+(1-dicodelta)*uncert)
}
