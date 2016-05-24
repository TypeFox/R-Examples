# {{{ header
#' Summary method for prodlim objects.
#' 
#' Summarizing the result of the product limit method in life-table format.
#' Calculates the number of subjects at risk and counts events and censored
#' observations at specified times or in specified time intervals.
#' 
#' For cluster-correlated data the number of clusters at-risk are are also
#' given. Confidence intervals are displayed when they are part of the fitted
#' object.
#' 
#' @param object An object with class `prodlim' derived with
#' \code{\link{prodlim}}
#' @param times Vector of times at which to return the estimated
#' probabilities.
#' @param newdata A data frame with the same variable names as those
#' that appear on the right hand side of the 'prodlim' formula.
#' Defaults to \code{object$X}.
#' @param max.tables Integer. If \code{newdata} is not given the value
#' of \code{max.tables} decides about the maximal number of tables to
#' be shown.  Defaults to 20.
#' @param surv Logical. If FALSE report event probabilities instead of
#' survival probabilities. Only available for
#' \code{object$model=="survival"}.
#' @param cause The cause for predicting the cause-specific cumulative
#' incidence function in competing risk models.
#' @param intervals Logical. If TRUE count events and censored in
#' intervals between the values of \code{times}.
#' @param percent Logical. If TRUE all estimated values are multiplied
#' by 100 and thus interpretable on a percent scale.
#' @param showTime If \code{TRUE} evaluation times are put into a
#' column of the output table, otherwise evaluation times are shown as
#' rownames.
#' @param asMatrix Control the output format when there are multiple
#' life tables, either because of covariate strata or competing causes
#' or both.  If not missing and not FALSE, reduce multiple life tables
#' into a matrix with new columns \code{X} for covariate strata and
#' \code{Event} for competing risks.
#' @param ... Further arguments that are passed to the print
#' function.
#' @return A data.frame with the relevant information.
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{prodlim}}, \code{\link{summary.Hist}}
#' 
#' @keywords survival
##' @examples
##' 
##' library(lava)
##' set.seed(17)
##' m <- survModel()
##' distribution(m,~age) <- uniform.lvm(30,80)
##' distribution(m,~sex) <- binomial.lvm()
##' m <- categorical(m,~z,K=3)
##' regression(m,eventtime~age) <- 0.01
##' regression(m,eventtime~sex) <- -0.4
##' d <- sim(m,50)
##' d$sex <- factor(d$sex,levels=c(0,1),labels=c("female","male"))
##' d$Z <- factor(d$z,levels=c(1,0,2),labels=c("B","A","C"))
##' 
##' # Univariate Kaplan-Meier
##' # -----------------------------------------------------------------------------------------
##' fit0 <- prodlim(Hist(time,event)~1,data=d)
##' summary(fit0)
##' 
##' ## show survival probabilities as percentage and
##' ## count number of events within intervals of a
##' ## given time-grid:
##' summary(fit0,times=c(1,5,10,12),percent=TRUE,intervals=TRUE)
##' 
##' ## the result of summary has a print function
##' ## which passes ... to print and print.listof
##' sx <- summary(fit0,times=c(1,5,10,12),percent=TRUE,intervals=TRUE)
##' print(sx,digits=3)
##' 
##' ## show cumulative incidences (1-survival)
##' summary(fit0,times=c(1,5,10,12),surv=FALSE,percent=TRUE,intervals=TRUE)
##' 
##' # Stratified Kaplan-Meier
##' # -----------------------------------------------------------------------------------------
##' 
##' fit1 <- prodlim(Hist(time,event)~sex,data=d)
##' print(summary(fit1,times=c(1,5,10),intervals=TRUE,percent=TRUE),digits=3)
##' 
##' summary(fit1,times=c(1,5,10),asMatrix=TRUE,intervals=TRUE,percent=TRUE)
##' 
##' fit2 <- prodlim(Hist(time,event)~Z,data=d)
##' print(summary(fit2,times=c(1,5,10),intervals=TRUE,percent=TRUE),digits=3)
##' 
##' ## Continuous strata (Beran estimator)
##' # -----------------------------------------------------------------------------------------
##' fit3 <- prodlim(Hist(time,event)~age,data=d)
##' print(summary(fit3,
##'               times=c(1,5,10),
##'               newdata=data.frame(age=c(20,50,70)),
##'               intervals=TRUE,
##'               percent=TRUE),digits=3)
##' 
##' ## stratified Beran estimator
##' # -----------------------------------------------------------------------------------------
##' fit4 <- prodlim(Hist(time,event)~age+sex,data=d)
##' print(summary(fit4,
##'               times=c(1,5,10),
##'               newdata=data.frame(age=c(20,50,70),sex=c("female","male","male")),
##'               intervals=TRUE,
##'               percent=TRUE),digits=3)
##' 
##' print(summary(fit4,
##'               times=c(1,5,10),
##'               newdata=data.frame(age=c(20,50,70),sex=c("female","male","male")),
##'               intervals=TRUE,collapse=TRUE,
##'               percent=TRUE),digits=3)
##' 
##' ## assess results from summary
##' x <- summary(fit4,times=10,newdata=expand.grid(age=c(60,40,50),sex=c("male","female")))
##' cbind(names(x$table),do.call("rbind",lapply(x$table,round,2)))
##' 
##' x <- summary(fit4,times=10,newdata=expand.grid(age=c(60,40,50),sex=c("male","female")))
##' 
##' ## Competing risks: Aalen-Johansen
##' # -----------------------------------------------------------------------------------------
##' d <- SimCompRisk(30)
##' crfit <- prodlim(Hist(time,event)~X1,data=d)
##' summary(crfit,times=c(1,2,5))
##' summary(crfit,times=c(1,2,5),cause=1,intervals=TRUE)
##' summary(crfit,times=c(1,2,5),cause=1,asMatrix=TRUE)
##' summary(crfit,times=c(1,2,5),cause=1:2,asMatrix=TRUE)
##' 
##' 
##' # extract the actual tables from the summary 
##' sumfit <- summary(crfit,times=c(1,2,5),print=FALSE)
##' sumfit$table[[1]] # cause 1
##' sumfit$table[[2]] # cause 2
##' 
##' 
##' # '
#' @export 
summary.prodlim <- function(object,
                            times,
                            newdata,
                            max.tables=20,
                            surv=TRUE,
                            cause,
                            intervals=FALSE,
                            percent=FALSE,
                            showTime=TRUE,
                            asMatrix=FALSE,
                            ...) {
    # }}}
    # {{{  classify the situation
    cens.type <- object$cens.type         # uncensored, right or interval censored
    model <- object$model                 # survival, competing risks or multi-state
    ## cluster <- object$clustervar          # clustered data?
    cotype <- object$covariate.type       # no, discrete, continuous or both
    # }}}
    # {{{  times
    jump.times <- object$time
    if (missing(times) && (length(times <- jump.times) > 50)) 
        times <- quantile(sort(unique(jump.times)))
    times <- sort(unique(times))
    if (any(times>max(jump.times)))
        warning(call.=TRUE,
                immediate.=TRUE,
                paste("\n","Time(s) ",paste(times[times>max(jump.times)],collapse=", "),
                      " are beyond the maximal follow-up time ",max(jump.times),"\n"))
    ntimes <- length(times)
    # }}}
    # {{{ interval-censored
    if (cens.type=="intervalCensored"){
        ltab <- data.frame(time=paste("(",paste(signif(object$time[1,],2),
                               signif(object$time[2,],2),
                               sep="-"),"]",sep=""),
                           n.risk=signif(object$n.risk,2),
                           n.event=signif(object$n.event,2),
                           ##    n.lost=object$n.lost,
                           surv=object$surv)
    }
    else{
        # }}}
        # {{{ with covariates
        if (cotype>1){
            if (missing(newdata) || length(newdata)==0){
                X <- object$X
                if (NROW(X)>max.tables){
                    warning(call.=TRUE,immediate.=TRUE,paste("\nLife tables are available for",
                                           NROW(X),
                                           "different covariate constellations.\n",
                                           "Shown are the table corresponding to the first row in object$X,",
                                           "corresponding to the middle row (median of the number of rows in object$X) ",
                                           "and corresponding to the last row in object$X ...\n",
                                           "to see more tables use arguments `newdata' and `max.tables'\n"))
                    X <- X[c(1,round(median(1:NROW(X))),NROW(X)),,drop=FALSE]
                }
            } else{
                  X <- unique.data.frame(newdata)
                  if (NROW(X) < NROW(newdata))
                      warning("Returned is only one summary for each unique value in newdata.")
              }
        } else {
              X <- NULL
          }
        if (model=="survival") {
            stats <- list(c("surv",1),c("se.surv",0))
            if (!is.null(object$conf.int))
                stats <- c(stats,list(c("lower",0),c("upper",1)))
            if (surv==FALSE){
                object$cuminc <- 1-object$surv
                object$se.cuminc <- object$se.surv
                cuminc.upper <- 1-object$lower
                cuminc.lower <- 1-object$upper
                object$lower <- cuminc.lower
                object$upper <- cuminc.upper
                stats <- list(c("cuminc",0),c("se.cuminc",0))
                if (!is.null(object$conf.int))
                    stats <- c(stats,list(c("lower",0),c("upper",1)))
            }
        }
        if (model=="competing.risks"){
            stats <- list(c("cuminc",0),c("se.cuminc",0))
            if (!is.null(object$conf.int))
                stats <- c(stats,list(c("lower",0),c("upper",0)))
            if (!missing(cause)){
                cause <- checkCauses(cause=cause,object=object)
            } else{ ## show all causes
                  cause <- attr(object$model.response,"states")
              }
            ltab <- lifeTab(object=object,
                            times=times,
                            cause=cause,
                            newdata=X,
                            stats=stats,
                            intervals=intervals,
                            percent=percent,
                            showTime=showTime)
            Found <- match(cause,names(ltab),nomatch=0)
            if (all(Found)>0) {
                ltab <- ltab[Found]
            }
            else stop(paste("\nCannot find cause: ",cause,".\nFitted were causes: ",paste(names(ltab),collapse=", "),sep=""))
        }else{
             ltab <- lifeTab(object=object,
                             times=times,
                             newdata=X,
                             stats=stats,
                             intervals=intervals,
                             percent=percent,
                             showTime=showTime)
         }
    }
    # }}}
    # {{{ output
    if (asMatrix!=FALSE) asMatrix <- TRUE
    if (model=="competing.risks"){
        ## out <- list(table=ltab,cause=cause)
        if (asMatrix)
            if (cotype>1)
                ltab <- List2Matrix(ltab,depth=2,names=c("Event","X"))
            else
                ltab <- List2Matrix(ltab,depth=1,names=c("Event"))
        
    }else{
         if(cotype>1 && asMatrix)
             ltab <- List2Matrix(ltab,depth=1,names="X")
     }
    out <- list(table=ltab,model=model,cotype=cotype,asMatrix=asMatrix,percent=percent)
    if (model=="competing.risks"){
        out <- c(out,list(cause=cause))
    }
    class(out) <- "summary.prodlim"
    out
    # }}}
}
