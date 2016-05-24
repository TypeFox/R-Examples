#' Survival model for simulation
#' 
#' Create a survival model to simulate a right censored event time data without
#' covariates
#' 
#' This function requires the \code{lava} package.
#' 
#' @return A structural equation model initialized with three variables: the
#' latent event time, the latent right censored time, and the observed
#' right censored event time.
#' @author Thomas A. Gerds <tag@@biostat.ku.dk>
#' @export 
survModel <- function(){
    ## require(lava)
    sm <- lava::lvm(~eventtime+censtime)
    lava::distribution(sm,"eventtime") <- lava::coxWeibull.lvm(scale=1/100)
    lava::distribution(sm,"censtime") <- lava::coxWeibull.lvm(scale=1/100)
    sm <- lava::eventTime(sm,time~min(eventtime=1,censtime=0),"event")
    sm
}
