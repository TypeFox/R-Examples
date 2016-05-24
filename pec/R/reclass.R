##' Retrospective table of risks predicted by two different methods, models, algorithms
##'
##' All risks are multiplied by 100 before 
##' @title Retrospective risk reclassification table
##' @param object Either a 
##' list with two elements. Each element should either
##' be a vector with probabilities, or an object for which
##' \code{predictSurvProb} or \code{predictEventProb} can extract predicted risk based on data.
##' @param reference Reference prediction model.
##' @param formula A survival formula as obtained either with
##' \code{prodlim::Hist} or \code{survival::Surv} which defines the
##' response in the \code{data}.
##' @param data Used to extract the response from the data and passed
##' on to \code{predictEventProb} to extract predicted event
##' probabilities.
##' @param time Time interest for prediction.
##' @param cause For competing risk models the cause of
##' interest. Defaults to all available causes.
##' @param cuts Risk quantiles to group risks.
##' @param digits Number of digits to show for the predicted risks
##' @return reclassification tables: overall table and one conditional table for each cause and for subjects event free at time interest.
##' @seealso predictStatusProb
##' @examples
##' \dontrun{
##' library(survival)
#' set.seed(40)
#' d <- prodlim::SimSurv(400)
#' nd <- prodlim::SimSurv(400)
#' Models <- list("Cox.X2"=coxph(Surv(time,status)~X2,data=d),
#'                "Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=d))
#' rc <- reclass(Models,formula=Surv(time,status)~1,data=nd,time=5)
#' print(rc)
#' plot(rc)
#'
#' set.seed(40)
#' library(riskRegression)
#' library(prodlim)
#' dcr <- prodlim::SimCompRisk(400)
#' ndcr <- prodlim::SimCompRisk(400)
#' crPred5 <- list("X2"=predictEventProb(CSC(Hist(time,event)~X2,data=dcr),newdata=ndcr,times=5),
#'                 "X1+X2"=predictEventProb(CSC(Hist(time,event)~X1+X2,data=dcr),newdata=ndcr,times=5))
#' rc <- reclass(crPred5,Hist(time,event)~1,data=ndcr,time=3)
#' print(rc)
#' 
#' reclass(crPred5,Hist(time,event)~1,data=ndcr,time=5,cuts=100*c(0,0.05,0.1,0.2,1))
#'}
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
reclass <- function(object,
                    reference,
                    formula,
                    data,
                    time,
                    cause,
                    cuts=seq(0,100,25),digits=2){
    if (missing(reference)){
        stopifnot(length(object)==2)
    } else{
          object <- list(object,reference)
      }
    NC <- length(cuts)
    NR <- NC-1 ## dimension of reclassification tables is NR x NR
    # {{{ response
    ## histformula <- formula
    ## if (histformula[[2]][[1]]==as.name("Surv")){
    ## histformula <- update(histformula,paste("prodlim::Hist","~."))
    ## histformula[[2]][[1]] <- as.name("prodlim::Hist")
    ## }
    ## print(histformula)
    ## m <- model.frame(histformula,data,na.action=na.fail)
    m <- model.frame(formula,data,na.action=na.omit)
    response <- model.response(m)
    if (match("Surv",class(response),nomatch=0)!=0){
        attr(response,"model") <- "survival"
        attr(response,"cens.type") <- "rightCensored"
        model.type <- "survival"
    }
    model.type <- attr(response,"model")
    if (model.type=="competing.risks"){
        predictHandlerFun <- "predictEventProb"
        availableCauses <- attr(response,"states")
        ncauses <- length(availableCauses)
        if (missing(cause))
            cause <- availableCauses[[1]]
        else
            if (match(cause, availableCauses,nomatch=FALSE)==0)
                stop("Cause ",cause," is not among the available causes: ",paste(availableCauses,collapse=", "))
    }
    else{
        predictHandlerFun <- "predictSurvProb"
    }
    # }}}
    ## for competing risks find the cause of interest.
    cutP <- function(P,cuts){
        if (min(P)<min(cuts))
            stop("Smallest predicted risk is smaller than first cut.")
        if (max(P)>max(cuts))
            stop("Largest predicted risk is larger than last cut.")
        cut(P,cuts,
            include.lowest=TRUE,
            labels=paste(paste(cuts[-NC],cuts[-1],sep="-"),"%",sep=""))
    }
    getPredictions <- function(x){
        if (any(is.na(x))) stop("Missing values in object.")
        P <- switch(class(x)[[1]],
                    "factor"={x},
                    "numeric"={
                        if (all(x<1)){
                            warning("Assumed that predictions are given on the scale [0,1] and multiplied by 100.")
                            x*100
                        } else{
                              x
                          }
                    },
                    {if (predictHandlerFun=="predictEventProb"){
                         P <- 100*do.call(predictHandlerFun,list(x,newdata=data,times=time,cause=cause))
                     } else {
                           P <- 100*do.call(predictHandlerFun,list(x,newdata=data,times=time))
                       }
                     P})
    }
    predrisk <- lapply(object,getPredictions)
    names(predrisk) <- names(object)
    predriskCut <- lapply(predrisk,function(P){if (is.factor(P)) P else cutP(P,cuts)})
    ## overall reclassification table
    retab <- table(predriskCut[[1]],predriskCut[[2]])
    ## reclassification frequencies conditional on outcome
    edat <- data.frame(cbind(do.call("cbind",predriskCut),response))
    edat$event[edat$status==0] <- 0
    N <- NROW(edat)
    names(edat)[1:2] <- c("P1","P2")
    cells <- split(edat,list(factor(edat$P1,levels=1:NR),factor(edat$P2,levels=1:NR)))
    all.comb <- apply(expand.grid(1:(NR),1:(NR)),1,paste,collapse=".")
    nn <- names(object)
    if (!is.null(nn) & length(nn)==2){
        names(dimnames(retab)) <- nn
    }
    ## --------------------------------------------------------------------------------------
    ## Apply Bayes' theorem to calculate expected reclassification probabilities
    ##                      conditional on outcome
    ## --------------------------------------------------------------------------------------
    if (predictHandlerFun=="predictEventProb"){
        ## --------------------------------------------------------------------------------------
        ## Competing risk
        ##
        ## P(X=x|T<=t, cause=j) = P(X=x,T<=t,cause=j)        / P(T<=t,cause=j)
        ##                      = P(T<=t,cause=j|X=x) P(X=x) / P(T<=t,cause=j)
        ##                      = cuminc.x H.x / cuminc
        ## 
        ##           P(X=x|T>t) = P(X=x,T>t) / P(T>t)
        ##                      = P(T>t|X=x) P(X=x) / P(T>t)
        ##                      = efreesurv.x H.x / efreesurv
        ## --------------------------------------------------------------------------------------
        eformula <- Hist(time,event)~1
        Hx <- unlist(lapply(cells,NROW))/N
        cuminc.x <- do.call("rbind",
                            lapply(names(cells),
                                   function(cc){
                                       x <- cells[[cc]]
                                       if (NROW(x)>0){
                                           ## warn if too short followup
                                           if (all(x$time<time)) {
                                               warning(call.=FALSE,paste0("pec::reclass: Cell row ",
                                                           sub("\\."," column ",cc),
                                                           " no subject was followed until time ",
                                                           time,
                                                           ". Result is NA (not available)."))
                                               rep(NA,length(availableCauses))
                                           } else{
                                                 fit.x <- prodlim::prodlim(eformula,data=x)
                                                 fitted.causes <- attr(fit.x$model.response,"states")
                                                 nstates <- length(fitted.causes)
                                                 sapply(availableCauses,function(j){
                                                            ## it may happen that cause j
                                                            ## does not occur in this cell
                                                            if (sum(x$event==j)>0){
                                                                ## check if there is more than one cause
                                                                if (nstates<length(availableCauses)){
                                                                    if (nstates==1){
                                                                        ## only one cause
                                                                        predict(fit.x,times=time,type="cuminc")
                                                                    } else{
                                                                          ## competing causes but less than all causes
                                                                          ## need to change the value of cause
                                                                          xj.cause <- match(j,fitted.causes,nomatch=0)
                                                                          if (xj.cause==0)
                                                                              stop(paste0("Cause ",j,"does not appear in fit. Fitted are causes: ",fitted.causes))
                                                                          else{
                                                                              predict(fit.x,times=time,cause=xj.cause,type="cuminc")
                                                                          }
                                                                      }
                                                                }else{
                                                                     predict(fit.x,times=time,cause=j,type="cuminc")
                                                                 }
                                                            } else {
                                                                  ## warn if no event of type j
                                                                  jstring <- j
                                                                  if (as.character(j)%in%availableCauses[[j]])
                                                                      jstring <- paste0(j," (",availableCauses[[j]],")")
                                                                  warning(call.=FALSE,paste0("pec::reclass: Cell row ",
                                                                              sub("\\."," column ",cc),
                                                                              " no event of type ",
                                                                              jstring,". Result is 0."))
                                                                  return(0)
                                                              }
                                                        })
                                             }
                                       } else{
                                             ## empty cell
                                             rep(0,length(availableCauses))
                                         }}))
        fit <- prodlim::prodlim(eformula,data=edat)
        cuminc <- unlist(lapply(availableCauses,function(j){predict(fit,times=time,cause=j,type="cuminc")}))
        Px <- apply(cuminc.x * Hx,1,function(p){p/ cuminc})
        ## rownames(Px) <- paste("cause",availableCauses,sep=":")
        efreesurv <- 1-sum(cuminc)
        efreesurv.x <- 1-rowSums(cuminc.x,na.rm=TRUE)
        Px <- rbind(Px,"eventfree"=efreesurv.x * Hx / efreesurv)
        event.retab <- lapply(1:NROW(Px),function(xx){
                                        matrix(Px[xx,],ncol=NR)
                                        matrix(Px[xx,],ncol=NR,dimnames=dimnames(retab))
                                    })
        names(event.retab) <- c(paste("Event:",availableCauses),"Event-free")
    } else{
          ## --------------------------------------------------------------------------------------
          ## Survival
          ##
          ## P(X=x|T<=t) = P(X=x,T<=t) /P(T<=t)
          ##             = P(T<=t|X=x) P(X=x) /P(T<=t)
          ##             = cuminc.x * Hx / cuminc
          ##
          ## P(X=x|T>t)  = P(X=x,T>t) /P(T>t)
          ##             = surv.x * Hx / surv
          ## 
          ## --------------------------------------------------------------------------------------
          eformula <- Hist(time,status)~1
          Hx <- unlist(lapply(cells,NROW))/N
          cuminc <- predict(prodlim::prodlim(eformula,data=edat),times=time,type="cuminc")
          cuminc.x <- sapply(cells,function(x){
                                 if (NROW(x)>0){
                                     ## warn if too short followup
                                     if (all(x$time<time)) {
                                         warning(call.=FALSE,paste0("pec::reclass: ",
                                                     ## sub("\\."," column ",cc),
                                                     " no subject was followed until time ",
                                                     time,
                                                     ". Result is NA (not available)."))
                                         NA
                                     }else{
                                          predict(prodlim::prodlim(eformula,data=x),times=time)
                                      }
                                 } else {
                                       ## empty cell
                                       0
                                   }
                             })
          surv <- 1-cuminc
          surv.x <- 1-cuminc.x
          dimnames=dimnames(retab)
          event.retab <- list("event"=matrix(cuminc.x*Hx/surv,ncol=NR),
                              "eventfree"=matrix(surv.x * Hx / surv,ncol=NR,dimnames=dimnames(retab)))
      }
    out <- list(time=time,
                predictedRisk=predrisk,
                reclassification=retab,
                event.reclassification=event.retab,
                cuts=cuts,
                model=model.type)
    class(out) <- "riskReclassification"
    out
}




