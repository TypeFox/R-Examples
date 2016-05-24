##' Simulate right censored survival data with two covariates X1 and X2, both have effect exp(1) on the hazard of the unobserved event time.
##'
##' This function calls \code{survModel}, then adds  covariates and finally calls \code{sim.lvm}.
##' @title Simulate survival data
##' @param N sample size
##' @param ... do nothing
##' @return data.frame with simulated data
##' @references Bender, Augustin & Blettner. Generating survival times to simulate Cox proportional hazards models. Statistics in Medicine, 24: 1713-1723, 2005.
##' @author Thomas Alexander Gerds
##' @examples
##' 
##' SimSurv(10)
##'
##' @export
SimSurv <- function(N, ...){
    m <- survModel()
    regression(m,from="X1",to="eventtime") <- 1
    regression(m,from="X2",to="eventtime") <- 1
    distribution(m,"X1") <- binomial.lvm()
    m <- eventTime(m,time~min(eventtime=1,censtime=0),"status")
    sim(m,N)
}

## SimSurvInternalIntervalCensored <- function(N,
                                            ## unit,
                                            ## lateness,
                                            ## compliance,
                                            ## withdraw.time,
                                            ## event.time){
  ## Intervals <- do.call("rbind",lapply(1:N,function(i){
    ## schedule <- seq(0,withdraw.time[i],unit)
    ## M <- length(schedule)
    ## g <- c(0,rep(unit,M))
    ## # introduce normal variation of the visit times
    ## g <- g+c(abs(rnorm(1,0,lateness)),rnorm(M,0,lateness))
    ## grid <- c(0,cumsum(g))
    ## # remove visits after the end of follow-up time
    ## grid <- grid[grid<withdraw.time[i]]
    ## # remove intermediate visits
    ## if (compliance<1){
      ## stopifnot(compliance>0)
      ## missed <- rbinom(length(grid),1,compliance)==0
      ## grid <- grid[missed==FALSE]
    ## }
    ## if (length(grid)==0){
      ## L <- 0
      ## R <- Inf
    ## }
    ## else{
      ## posTime <- sindex(jump.times=grid,
                        ## eval.times=event.time[i])
      ## L <- grid[posTime]
      ## R <- grid[posTime+1]
      ## if (is.na(R)){
        ## R <- Inf
      ## }
    ## }
    ## c(L=L,R=R)
  ## }))
  ## out <- data.frame(Intervals)
  ## out
## }
# }}}
# {{{ find.baseline
## find.baseline <- function(x=.5,
                          ## setting,
                          ## verbose=FALSE){
  ## N <- setting$N
  ## f <- function(y){
    ## setting$cens.baseline <- y
    ## ncens <- sum(do.call("SimSurv",replace(setting,"verbose",verbose))$status==0)
    ## x-ncens/N
  ## }
  ## base.cens <- uniroot(f,c(exp(-50),1000000),tol=.0000001,maxiter=100)$root
  ## new.setting <- setting
  ## new.setting$cens.baseline <- base.cens
  ## do.call("SimSurv",replace(new.setting,"verbose",TRUE))
  ## new.setting
## }
# }}}
# {{{quantile.SimSurv
## quantile.SimSurv <- function(x,B=10,na.rm=FALSE,probs=.9){
  ## callx <- attr(x,"call")
  ## nix <- do.call("rbind",lapply(1:B,function(b){
    ## quantile(eval(callx)$time,probs)
  ## }))
  ## nix <- colMeans(nix)
  ## nix
## }
# }}}
