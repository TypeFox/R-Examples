##' Extract event history data and design matrix including specials from call
##'
##' Obtain a list with the data used for event history regression analysis. This
##' function cannot be used directly on the user level but inside a function
##' to prepare data for survival analysis. 
##' @title Event history frame
##' @param formula Formula whose left hand side specifies the event
##' history, i.e., either via Surv() or Hist().
##' @param data Data frame in which the formula is interpreted
##' @param unspecialsDesign Passed as is to
##' \code{\link{model.design}}.
##' @param specials Character vector of special function names.
##' Usually the body of the special functions is function(x)x but
##' e.g., \code{\link{strata}} from the survival package does treat
##' the values
##' @param specialsFactor Passed as is to \code{\link{model.design}}.
##' @param specialsDesign Passed as is to \code{\link{model.design}}
##' @param stripSpecials Passed as \code{specials} to
##' \code{\link{strip.terms}}
##' @param stripArguments Passed as \code{arguments} to
##' \code{\link{strip.terms}}
##' @param stripAlias Passed as \code{alias.names} to
##' \code{\link{strip.terms}}
##' @param stripUnspecials Passed as \code{unspecials} to
##' \code{\link{strip.terms}}
##' @param dropIntercept Passed as is to \code{\link{model.design}}
##' @param check.formula If TRUE check if formula is a Surv or Hist
##' thing.
##' @param response If FALSE do not get response data (event.history).
##' @return A list which contains
##' - the event.history (see \code{\link{Hist}})
##' - the design matrix (see \code{\link{model.design}})
##' - one entry for each special (see \code{\link{model.design}})
##' @seealso model.frame model.design Hist
##' @examples
##' 
##' ## Here are some data with an event time and no competing risks
##' ## and two covariates X1 and X2.
##' ## Suppose we want to declare that variable X1 is treated differently
##' ## than variable X2. For example, X1 could be a cluster variable, or
##' ## X1 should have a proportional effect on the outcome.
##' dsurv <- data.frame(time=1:7,
##'                     status=c(0,1,1,0,0,0,1),
##'                     X2=c(2.24,3.22,9.59,4.4,3.54,6.81,5.05),
##'                     X3=c(1,1,1,1,0,0,1),
##'                     X4=c(44.69,37.41,68.54,38.85,35.9,27.02,41.84),
##'                     X1=factor(c("a","b","a","c","c","a","b"),
##'                         levels=c("c","a","b")))
##' ## We pass a formula and the data
##' e <- EventHistory.frame(Hist(time,status)~prop(X1)+X2+cluster(X3)+X4,
##'                         data=dsurv,
##'                         specials=c("prop","cluster"),
##'                         stripSpecials=c("prop","cluster"))
##' names(e)
##' ## The first element is the event.history which is result of the left hand
##' ## side of the formula:
##' e$event.history
##' ## same as
##' with(dsurv,Hist(time,status))
##' ## to see the structure do 
##' colnames(e$event.history)
##' unclass(e$event.history)
##' ## in case of competing risks there will be an additional column called event,
##' ## see help(Hist) for more details
##' 
##' ## The other elements are the design, i.e., model.matrix for the non-special covariates
##' e$design
##' ## and a data.frame for the special covariates
##' e$prop
##' ## The special covariates can be returned as a model.matrix
##' e2 <- EventHistory.frame(Hist(time,status)~prop(X1)+X2+cluster(X3)+X4,
##'                          data=dsurv,
##'                          specials=c("prop","cluster"),
##'                          stripSpecials=c("prop","cluster"),
##'                          specialsDesign=TRUE)
##' e2$prop
##' ## and the non-special covariates can be returned as a data.frame
##' e3 <- EventHistory.frame(Hist(time,status)~prop(X1)+X2+cluster(X3)+X4,
##'                          data=dsurv,
##'                          specials=c("prop","cluster"),
##'                          stripSpecials=c("prop","cluster"),
##'                          specialsDesign=TRUE,
##'                          unspecialsDesign=FALSE)
##' e3$design
##' 
##' ## the general idea is that the function is used to parse the combination of
##' ## formula and data inside another function. Here is an example with
##' ## competing risks
##' SampleRegression <- function(formula,data=parent.frame()){
##'     thecall <- match.call()
##'     ehf <- EventHistory.frame(formula=formula,
##'                               data=data,
##'                               stripSpecials=c("prop","cluster","timevar"),
##'                               specials=c("prop","timevar","cluster"))
##'     time <- ehf$event.history[,"time"]
##'     status <- ehf$event.history[,"status"]
##'     ## event as a factor
##'     if (attr(ehf$event.history,"model")=="competing.risks"){
##'         event <- ehf$event.history[,"event"]
##'         Event <- getEvent(ehf$event.history)
##'         list(response=data.frame(time,status,event,Event),X=ehf[-1])
##'     }
##'     else{ # no competing risks
##'         list(response=data.frame(time,status),X=ehf[-1])
##'     }
##' }
##' dsurv$outcome <- c("cause1","0","cause2","cause1","cause2","cause2","0")
##' SampleRegression(Hist(time,outcome)~prop(X1)+X2+cluster(X3)+X4,dsurv)
##' 
##' ## let's test if the parsing works
##' form1 <- Hist(time,outcome!="0")~prop(X1)+X2+cluster(X3)+X4
##' form2 <- Hist(time,outcome)~prop(X1)+cluster(X3)+X4
##' ff <- list(form1,form2)
##' lapply(ff,function(f){SampleRegression(f,dsurv)})
##' 
##' 
##' ## here is what the riskRegression package uses to
##' ## distinguish between covariates with
##' ## time-proportional effects and covariates with
##' ## time-varying effects:
##' \dontrun{
##' library(riskRegression)
##' data(Melanoma)
##' f <- Hist(time,status)~prop(thick)+strata(sex)+age+prop(ulcer,power=1)+timevar(invasion,test=1)
##' ## here the unspecial terms, i.e., the term age is treated as prop
##' ## also, strata is an alias for timvar
##' 
##' EHF <- prodlim::EventHistory.frame(formula,
##'                                    Melanoma[1:10],
##'                                    specials=c("timevar","strata","prop","const","tp"),
##'                                    stripSpecials=c("timevar","prop"),
##'                                    stripArguments=list("prop"=list("power"=0),
##'                                        "timevar"=list("test"=0)),
##'                                    stripAlias=list("timevar"=c("strata"),
##'                                        "prop"=c("tp","const")),
##'                                    stripUnspecials="prop",
##'                                    specialsDesign=TRUE,
##'                                    dropIntercept=TRUE)       
##' EHF$prop
##' EHF$timevar
##' }
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
EventHistory.frame <- function(formula,
                               data,
                               unspecialsDesign=TRUE,
                               specials,
                               specialsFactor=TRUE,
                               specialsDesign=FALSE,
                               stripSpecials=NULL,
                               stripArguments=NULL,
                               stripAlias=NULL,
                               stripUnspecials=NULL,
                               dropIntercept=TRUE,
                               check.formula=TRUE,
                               response=TRUE){
    # {{{  check if formula is a proper formula
    if (response && check.formula){
        formula.names <- try(all.names(formula),silent=TRUE)
        if (!(formula.names[1]=="~")
            ||
            (match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
            stop("Invalid specification of formula. Perhaps forgotten right hand side?\nNote that any subsetting, ie data$var or data[,\"var\"], is invalid for this function.")}
        else{
            if (!(any(match(c("survival::Surv","Surv","prodlim::Hist","Hist"),
                            formula.names,nomatch=0))))
                stop("formula is NOT a proper survival formula,\nwhich must have a `Surv' or `Hist' object as response.")
        }
    }
    # }}}
    # {{{call model.frame
    ## data argument is used to resolve '.' see help(terms.formula)
    Terms <- terms(x=formula,specials=specials,data=data)
    if (!is.null(stripSpecials)){
        ## Terms <- terms(x=formula, specials=specials)
        if (length(attr(Terms,"term.labels"))>0)
            Terms <- strip.terms(Terms,
                                 specials=stripSpecials,
                                 arguments=stripArguments,
                                 alias.names=stripAlias,
                                 unspecials=stripUnspecials)
    }
    # }}}
    # {{{ get all variables and remove missing values
    ## use the stripped formula because, otherwise
    ## it may be hard to know what variables are, e.g.,
    ## FGR uses cov2(var,tf=qfun) where qfun is a function 
    mm <- na.omit(get_all_vars(formula(Terms),data))
    if (NROW(mm) == 0) stop("No (non-missing) observations")
    # }}}

    # {{{ extract response
    if (response==TRUE && attr(Terms,"response")!=0){
        event.history <- model.response(model.frame(update(formula,".~1"),
                                                    data=mm))
        # }}}
        # {{{ Fix for those who use `Surv' instead of `Hist' 
        if (match("Surv",class(event.history),nomatch=0)!=0){
            attr(event.history,"model") <- "survival"
            attr(event.history,"cens.type") <- "rightCensored"
            attr(event.history,"entry.type") <- ifelse(ncol(event.history)==2,"","leftTruncated")
            if (attr(event.history,"entry.type")=="leftTruncated")
                colnames(event.history) <- c("entry","time","status")
        }
        # }}}
    }else event.history <- NULL
    # {{{ design
    design <- model.design(Terms,
                           data=mm,
                           maxOrder=1,
                           dropIntercept=dropIntercept,
                           unspecialsDesign=unspecialsDesign,
                           specialsFactor=specialsFactor,
                           specialsDesign=specialsDesign)
    # }}}
    out <- c(list(event.history=event.history),
             design[sapply(design,length)>0])
    attr(out,"Terms") <- Terms
    attr(out,"na.action") <- attr(mm,"na.action")
    class(out) <- "EventHistory.frame"
    out
}
##' @export 
as.data.frame.EventHistory.frame <- function(x,...){
    Y <- data.frame(unclass(x$event.history))
    X <- do.call("cbind",x[-1])
    cbind(Y,X)
}
