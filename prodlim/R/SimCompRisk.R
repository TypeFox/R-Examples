##' Simulate right censored competing risks data with two covariates X1 and X2. Both covariates have effect exp(1) on the hazards of event 1 and zero effect on the hazard of event 2.
##'
##' This function calls \code{crModel}, then adds covariates and finally calls \code{sim.lvm}.
##' @title Simulate competing risks data
##' @param N sample size
##' @param ... do nothing.
##' @return data.frame with simulated data
##' @author Thomas Alexander Gerds
##' @examples
##' 
##' SimCompRisk(10)
##'             
##' @export
SimCompRisk <- function(N, ...){
    ## require(lava)
    m <- crModel()
    regression(m,from="X1",to="eventtime1") <- 1
    regression(m,from="X2",to="eventtime1") <- 1
    distribution(m,"X1") <- binomial.lvm()
    out <- sim(m,N)
    ## for backward compatibility
    out$cause <- out$event
    out
}
 
