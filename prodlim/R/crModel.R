#' Competing risks model for simulation
#' 
#' Create a competing risks model with to causes to simulate a right censored event time data without
#' covariates
#' 
#' This function requires the \code{lava} package.
#' @title Competing risks model for simulation
#' @return A structural equation model initialized with four variables: the
#' latent event times of two causes, the latent right censored time, and the observed
#' right censored event time.
#' @author Thomas A. Gerds
#' @examples
#' library(lava)
#' m <- crModel()
#' d <- sim(m,6)
#' print(d)
#'
#' @export
crModel <- function(){
    # require(lava)
    crm <- lava::lvm()
    lava::distribution(crm,"eventtime1") <- lava::coxWeibull.lvm(scale=1/100)
    lava::distribution(crm,"eventtime2") <- lava::coxWeibull.lvm(scale=1/100)
    lava::distribution(crm,"censtime") <- lava::coxWeibull.lvm(scale=1/100)
    crm <- lava::eventTime(crm,time~min(eventtime1=1,eventtime2=2,censtime=0),"event")
    crm
}
