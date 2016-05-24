### stopTime.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Nov 28 2015 (10:07) 
## Version: 
## last-updated: Dec  4 2015 (06:57) 
##           By: Thomas Alexander Gerds
##     Update #: 23
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
##' All event times are stopped at a given time point and
##' corresponding events are censored
##'
##' @title Stop the time of an event history object
##' @param object Event history object as obtained with \code{Hist}
##' @param stop.time Time point at which to stop the event history object
##' @return Stopped event history object where all times are censored 
##'         at \code{stop.time}. All observations with times greater than \code{stop.time}
##'         are set to \code{stop.time} and the event status is set to \code{attr(object,"cens.code")}.
##'         A new column \code{"stop.time"} is equal to \code{1} for stopped observations
##'         and equal to \code{0} for the other observations.
##' @seealso Hist 
##' @examples
##'
##' set.seed(29)
##' d <- SimSurv(10)
##' h <- with(d,Hist(time,status))
##' h
##' stopTime(h,8)
##' stopTime(h,5)
##'
##' ## works also with Surv objects
##' library(survival)
##' s <- with(d,Surv(time,status))
##' stopTime(s,5)
##'
##' ## competing risks
##' set.seed(29)
##' dr <- SimCompRisk(10)
##' hr <- with(dr,Hist(time,event))
##' hr
##' stopTime(hr,8)
##' stopTime(hr,5)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
stopTime <- function(object,stop.time){
    if (missing(stop.time)) stop("Argument stop.time missing. Need a time point at which to stop the event history.")
    if (length(stop.time)>1) {
        warning("Argument stop.time is a vector. Proceed with the first element.")
        stop.time <- stop.time[[1]]
    }
    cc <- class(object)[[1]]
    stopifnot(cc%in% c("Hist","Surv"))
    if (cc=="Surv"){
        model <- "survival"
    }else{
         model <- attr(object,"model")
         if(!(model %in% c("survival","competing.risks")))
             stop(paste("Don't know (not yet) how to stop this type of model:",model))
     }
    stopped <- object[,"time"] >= stop.time
    sobject <- cbind(object,"stopped"=1*stopped)
    sobject[,"status"][stopped] <- 0
    if(model=="competing.risks") 
        sobject[,"event"][stopped] <- length(attr(object,"states"))+1
    sobject[,"time"][stopped] <- stop.time
    attr(sobject,"stop.time") <- stop.time
    attr(sobject,"class") <- attr(object,"class")
    if (cc=="Surv"){
        attr(sobject,"type") <- attr(object,"type")
    }
    attr(sobject,"states") <- attr(object,"states")
    attr(sobject,"model") <- attr(object,"model")
    attr(sobject,"cens.type") <- attr(object,"cens.type")
    attr(sobject,"cens.code") <- attr(object,"cens.code")
    attr(sobject,"entry.type") <- attr(object,"entry.type")
    sobject
}
#----------------------------------------------------------------------
### stopTime.R ends here
