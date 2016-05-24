#' @title Simulation routine for CUB models with shelter effect
#' @aliases simcubshe
#' @description Generate \eqn{n} pseudo-random numbers according to the CUB distribution
#'  with shelter effect and with given parameters.
#' @keywords distribution
#' @usage simcubshe(n, m, pai, csi, delta, shelter)
#' @export simcubshe
#' @import stats
#' @param n Number of simulated observations
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter
#' @param delta Shelter parameter
#' @param shelter Category corresponding to the shelter choice
#' @seealso \code{\link{probcubshe1}}, \code{\link{probcubshe2}}, \code{\link{probcubshe3}}
#' @examples
#' n<-300
#' m<-9
#' pai<-0.7
#' csi<-0.3
#' delta<-0.2
#' shelter<-3
#' simulation<-simcubshe(n, m, pai, csi, delta, shelter)
#' plot(table(simulation), xlab="Ordinal categories",ylab="Frequencies")


simcubshe <-
function(n,m,pai,csi,delta,shelter){
  dicopai<-runif(n)<pai
  dicodelta<-runif(n)<delta
  cub00<-dicopai*(1+rbinom(n,m-1,1-csi))+(1-dicopai)*sample(m,n,replace=TRUE)
  return((1-dicodelta)*cub00+dicodelta*shelter)
}
