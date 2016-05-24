##' product limit method
##' 
##' Nonparametric estimation in event history analysis. Featuring fast
##' algorithms and user friendly syntax adapted from the survival package.  The
##' product limit algorithm is used for right censored data; the
##' self-consistency algorithm for interval censored data.
##' 
##' 
##' The response of \code{formula} (ie the left hand side of the `~' operator)
##' specifies the model.
##' 
##' In two-state models -- the classical survival case -- the standard
##' Kaplan-Meier method is applied.  For this the response can be specified as a
##' \code{\link{Surv}} or as a \code{\link{Hist}} object. The \code{\link{Hist}}
##' function allows you to change the code for censored observations, e.g.
##' \code{Hist(time,status,cens.code="4")}.
##' 
##' Besides a slight gain of computing efficiency, there are some extensions
##' that are not included in the current version of the survival package:
##' 
##' (0) The Kaplan-Meier estimator for the censoring times \code{reverse=TRUE}
##' is correctly estimated when there are ties between event and censoring
##' times.
##' 
##' (1) A conditional version of the kernel smoothed Kaplan-Meier estimator for at most one
##' continuous predictors using nearest neighborhoods (Beran 1981,
##' Stute 1984, Akritas 1994).
##' 
##' (2) For cluster-correlated data the right hand side of \code{formula} may
##' identify a \code{\link{cluster}} variable. In that case Greenwood's variance
##' formula is replaced by the formula of Ying \& Wei (1994).
##' 
##' (3) Competing risk models can be specified via \code{\link{Hist}} response
##' objects in \code{formula}.
##' 
##' The Aalen-Johansen estimator is applied for estimating the cumulative
##' incidence functions for all causes.  The advantage over the function
##' \code{cuminc} of the cmprsk package are user-friendly model specification
##' via \code{\link{Hist}} and sophisticated print, summary, predict and plot
##' methods.
##' 
##' Under construction:
##' 
##' (U0) Interval censored event times specified via \code{\link{Hist}} are used
##' to find the nonparametric maximum likelihood estimate. Currently this works
##' only for two-state models and the results should match with those from the
##' package `Icens'.
##' 
##' (U1) Extensions to more complex multi-states models
##' 
##' (U2) The nonparametric maximum likelihood estimate for interval censored
##' observations of competing risks models.
##' 
##' @param formula A formula whose left hand side is a \code{Hist}
##' object. In some special cases it can also be a \code{Surv}
##' response object, see the details section. The right hand side is
##' as usual a linear combination of covariates which may contain at
##' most one continuous factor. Whether or not a covariate is
##' recognized as continuous or discrete depends on its class and on
##' the argument \code{discrete.level}. The right hand side may also
##' be used to specify clusters, see the details section.
##' @param data A data.frame in which all the variables of
##' \code{formula} can be interpreted.
##' @param subset Passed as argument \code{subset} to function
##' \code{subset} which applied to \code{data} before the formula is
##' processed.
##' @param na.action All lines in data with any missing values in the
##' variables of formula are removed.
##' @param reverse For right censored data, if reverse=TRUE then the
##' censoring distribution is estimated.
##' @param conf.int The level (between 0 and 1) for two-sided
##' pointwise confidence intervals. Defaults to 0.95. Remark: only
##' plain Wald-type confidence limits are available.
##' @param bandwidth Smoothing parameter for nearest neighborhoods
##' based on the values of a continuous covariate. See function
##' \code{neighborhood} for details.
##' @param caseweights Weights applied to the contribution of each
##' subject to change the number of events and the number at
##' risk. This can be used for bootstrap and survey analysis. Should
##' be a vector of the same length and the same order as \code{data}.
##' @param discrete.level Numeric covariates are treated as factors
##' when their number of unique values exceeds not
##' \code{discrete.level}. Otherwise the product limit method is
##' applied, in overlapping neighborhoods according to the bandwidth.
##' @param x logical value: if \code{TRUE}, the full covariate matrix
##' with is returned in component \code{model.matrix}.  The reduced
##' matrix contains unique rows of the full covariate matrix and is
##' always returned in component \code{X}.
##' @param method For interval censored data only.  If equal to
##' \code{"npmle"} (the default) use the usual Turnbull algorithm,
##' else the product limit version of the self-consistent estimate.
##' @param exact If TRUE the grid of time points used for estimation
##' includes all the L and R endpoints of the observed intervals.
##' @param maxiter For interval censored data only.  Maximal number of
##' iterations to obtain the nonparametric maximum likelihood
##' estimate.  Defaults to 1000.
##' @param grid For interval censored data only. When method=one.step
##' grid for one-step product limit estimate. Defaults to sorted list
##' of unique left and right endpoints of the observed intervals.
##' @param tol For interval censored data only. Numeric value whose
##' negative exponential is used as convergence criterion for finding
##' the nonparametric maximum likelihood estimate.  Defaults to 7
##' meaning exp(-7).
##' @param type In two state models either \code{"surv"} for the Kaplan-Meier estimate of the survival
##' function or \code{"cuminc"} for 1-Kaplan-Meier. Default is \code{"surv"} when \code{reverse==FALSE} and \code{"cuminc"} when \code{reverse==TRUE}.
##' In competing risks models it has to be \code{"cuminc"}
##' Aalen-Johansen estimate of the cumulative incidence function.
##' @return Object of class "prodlim". See \code{\link{print.prodlim}}, \code{\link{predict.prodlim}}, predict,
##' \code{\link{summary.prodlim}}, \code{\link{plot.prodlim}}.
##' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
##' @seealso \code{\link{predictSurv}}, \code{\link{predictSurvIndividual}},
##' \code{\link{predictCuminc}}, \code{\link{Hist}}, \code{\link{neighborhood}},
##' \code{\link{Surv}}, \code{\link{survfit}}, \code{\link{strata}},
##' @references Andersen, Borgan, Gill, Keiding (1993) Springer `Statistical
##' Models Based on Counting Processes'
##' 
##' Akritas (1994) The Annals of Statistics 22, 1299-1327 Nearest neighbor
##' estimation of a bivariate distribution under random censoring.
##' 
##' R Beran (1981) http://anson.ucdavis.edu/~beran/paper.html `Nonparametric
##' regression with randomly censored survival data'
##' 
##' Stute (1984) The Annals of Statistics 12, 917--926 `Asymptotic Normality of
##' Nearest Neighbor Regression Function Estimates'
##' 
##' Ying, Wei (1994) Journal of Multivariate Analysis 50, 17-29 The Kaplan-Meier
##' estimate for dependent failure time observations
##' @keywords survival nonparametric cluster
##' @examples
##' 
##' ##---------------------two-state survival model------------
##' dat <- SimSurv(30)
##' with(dat,plot(Hist(time,status)))
##' fit <- prodlim(Hist(time,status)~1,data=dat)
##' print(fit)
##' plot(fit)
##' summary(fit)
##' quantile(fit)
##'
##' ## Subset
##' fit1a <- prodlim(Hist(time,status)~1,data=dat,subset=dat$X1==1)
##' fit1b <- prodlim(Hist(time,status)~1,data=dat,subset=dat$X1==1 & dat$X2>0)
##' 
##' ## --------------------clustered data---------------------
##' library(survival)
##' cdat <- cbind(SimSurv(30),patnr=sample(1:5,size=30,replace=TRUE))
##' fit <- prodlim(Hist(time,status)~cluster(patnr),data=cdat)
##' print(fit)
##' plot(fit)
##' summary(fit)
##' 
##' 
##' ##-----------compare Kaplan-Meier to survival package---------
##' 
##' dat <- SimSurv(30)
##' pfit <- prodlim(Surv(time,status)~1,data=dat)
##' pfit <- prodlim(Hist(time,status)~1,data=dat) ## same thing
##' sfit <- survfit(Surv(time,status)~1,data=dat,conf.type="plain")
##' ##  same result for the survival distribution function 
##' all(round(pfit$surv,12)==round(sfit$surv,12))
##' summary(pfit,digits=3)
##' summary(sfit,times=quantile(unique(dat$time)))
##' 
##' ##-----------estimating the censoring survival function----------------
##' 
##' rdat <- data.frame(time=c(1,2,3,3,3,4,5,5,6,7),status=c(1,0,0,1,0,1,0,1,1,0))
##' rpfit <- prodlim(Hist(time,status)~1,data=rdat,reverse=TRUE)
##' rsfit <- survfit(Surv(time,1-status)~1,data=rdat,conf.type="plain")
##' ## When there are ties between times at which events are observed
##' ## times at which subjects are right censored, then the convention
##' ## is that events come first. This is not obeyed by the above call to survfit,
##' ## and hence only prodlim delivers the correct reverse Kaplan-Meier:
##' cbind("Wrong:"=rsfit$surv,"Correct:"=rpfit$surv)
##' 
##' ##-------------------stratified Kaplan-Meier---------------------
##' 
##' pfit.X2 <- prodlim(Surv(time,status)~X2,data=dat)
##' summary(pfit.X2)
##' summary(pfit.X2,intervals=TRUE)
##' plot(pfit.X2)
##' 
##' ##----------continuous covariate: Stone-Beran estimate------------
##' 
##' prodlim(Surv(time,status)~X1,data=dat)
##' 
##' ##-------------both discrete and continuous covariates------------
##' 
##' prodlim(Surv(time,status)~X2+X1,data=dat)
##' 
##' ##----------------------interval censored data----------------------
##' 
##' dat <- data.frame(L=1:10,R=c(2,3,12,8,9,10,7,12,12,12),status=c(1,1,0,1,1,1,1,0,0,0))
##' with(dat,Hist(time=list(L,R),event=status))
##' 
##' dat$event=1
##' npmle.fitml <- prodlim(Hist(time=list(L,R),event)~1,data=dat)
##' 
##' ##-------------competing risks-------------------
##' 
##' CompRiskFrame <- data.frame(time=1:100,event=rbinom(100,2,.5),X=rbinom(100,1,.5))
##' crFit <- prodlim(Hist(time,event)~X,data=CompRiskFrame)
##' summary(crFit)
##' plot(crFit)
##' summary(crFit,cause=2)
##' plot(crFit,cause=2)
##' 
##' 
##' # Changing the cens.code:
##' dat <- data.frame(time=1:10,status=c(1,2,1,2,5,5,1,1,2,2))
##' fit <- prodlim(Hist(time,status)~1,data=dat)
##' print(fit$model.response)
##' fit <- prodlim(Hist(time,status,cens.code="2")~1,data=dat)
##' print(fit$model.response)
##' plot(fit)
##' plot(fit,cause="5")
##' 
##' 
##' ##------------delayed entry----------------------
##' 
##' ## left-truncated event times with competing risk endpoint 
##' 
##' dat <- data.frame(entry=c(7,3,11,12,11,2,1,7,15,17,3),time=10:20,status=c(1,0,2,2,0,0,1,2,0,2,0))
##' fitd <- prodlim(Hist(time=time,event=status,entry=entry)~1,data=dat)
##' summary(fitd)
##' plot(fitd)
##' 
#' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
"prodlim" <- function(formula,
                      data = parent.frame(),
                      subset,
                      na.action=NULL,
                      reverse=FALSE,
                      conf.int=.95,
                      bandwidth=NULL,
                      caseweights,
                      discrete.level=3,
                      x=TRUE,
                      # force.multistate=FALSE,
                      maxiter=1000,
                      grid,
                      tol=7,
                      method=c("npmle","one.step","impute.midpoint","impute.right"),
                      exact=TRUE,
                      type){
    # {{{  find the data
    call <- match.call()
    if (!missing(subset))
        data <- subset(data,subset=subset)
    EHF <- EventHistory.frame(formula=formula,
                              data=data,
                              unspecialsDesign=FALSE,
                              specials=c("Strata","strata","factor", "NN","cluster"),
                              stripSpecials=c("strata","cluster","NN"),
                              stripAlias=list("strata"=c("Strata","factor")),
                              stripArguments=list("strata"=NULL,"NN"=NULL,"cluster"=NULL),
                              specialsDesign=FALSE,
                              check.formula=TRUE)
    event.history <- EHF$event.history
    response <- EHF$event.history
    if (reverse==TRUE){ ## estimation of censoring distribution
        model.type <- 1
    }else{
         model.type <- match(attr(event.history,"model"),c("survival","competing.risks","multi.states"))
     }
    if (missing(type)) type <- switch(model.type,"survival"=ifelse(reverse,"cuminc","surv"),"cuminc")
    else {
        type <- tolower(type)
        stopifnot(match(type,c("surv","cuminc"),nomatch=0)!=0)
    }
    cens.type <- attr(response,"cens.type")
    #  if (force.multistate==TRUE) model.type <- 3
    # {{{ order according to event times
    if (cens.type!="intervalCensored"){
        event.time.order <- order(event.history[,"time"],-event.history[,"status"])
    }
    else{
        event.time.order <- order(event.history[,"L"],-event.history[,"status"])
    }
    # }}}
    # {{{  covariates

    covariates <- EHF[-1]
    ##  `factor' and 'Strata' are aliases for `strata'
    strata.pos <- match(c("strata","factor","Strata"),names(covariates),nomatch=0)
    if (sum(strata.pos)>0)
        strata <- do.call("cbind",covariates[strata.pos])
    else
        strata <- NULL
    ##  'NN'
    NN <- covariates[["NN"]]
    xlevels <- attr(strata,"levels")
    ## unspecified
    rest <- covariates$design
    xlevels <- c(attr(strata,"levels"),attr(rest,"levels"))
    if ((is.null(NN)+is.null(strata)+is.null(rest))==3){
        cotype <- 1
    } else{
        unspecified <- NULL
        if (!is.null(rest)){
            discrete.p <- sapply(colnames(rest),function(u){
                x <- rest[,u,drop=TRUE]
                !is.numeric(x) || !length(unique(x))>discrete.level
            })
            if (any(!discrete.p)){ ## continuous covariates
                NN <- if (is.null(NN))
                          rest[,!discrete.p,drop=FALSE]
                      else
                          cbind(NN,rest[,!discrete.p,drop=FALSE])
            }
            if (any(discrete.p)){ ## discrete covariates
                strata <- if (is.null(strata)){
                    rest[,discrete.p,drop=FALSE]
                } else{
                    cbind(strata,rest[,discrete.p,drop=FALSE])
                }
            }
        }
        if (NCOL(NN)>1) stop(paste("Currently we can not compute neighborhoods in",length(colnames(NN)),"continuous dimensions."))
        cotype <- 1 + (!is.null(strata))*1+(!is.null(NN))*2
    }
    ## use unique values as levels
    ## for variables that are not factors
    ## but treated as such
    if (any(found <- (match(colnames(strata),names(xlevels),nomatch=0)==0))){
        uniquelevels <- lapply(colnames(strata)[found],function(x){
            unique(strata[,x])
        })
        names(uniquelevels) <- colnames(strata)[found]
        xlevels <- c(xlevels,uniquelevels)
    }
    ## cotype
    # 1 : no covariates
    # 2 : only strata 
    # 3 : only continuous
    # 4 : strata AND continuous
    # }}}
    # {{{  disjunct strata (discrete covariates)
    if (cotype %in% c(2,4)){
        ## changed 09 Dec 2014 (16:57)-->
        ## S <- do.call("paste", c(data.frame(strata), sep = "\r"))
        S <- interaction(data.frame(strata), sep = ":",drop=TRUE)
        ## <-- changed 09 Dec 2014 (16:57) 
        NS <- length(unique(S))
        ## changed 09 Dec 2014 (16:57) -->
        Sfactor <- factor(S,levels=levels(S),labels=1:NS)
        ## <-- changed 09 Dec 2014 (16:57)
        if (cens.type!="intervalCensored"){
            sorted <- order(Sfactor, response[,"time"],-response[,"status"])
        }
        else{
            sorted <- order(Sfactor, response[,"L"],-response[,"status"])
        }
        Sfactor <- Sfactor[sorted]
    }
    else{
        sorted <- event.time.order
    }
  
    response <- response[sorted,] # sort each stratum
  
    # }}}
    # {{{  overlapping neighborhoods (continuous covariates)
  
    if (cotype %in% c(3,4)){
        Z <- NN[sorted,,drop=TRUE]
        if (cotype==3){
            nbh <- neighborhood(Z,bandwidth=bandwidth) 
            nbh.list <- list(nbh)
            bandwidth <- nbh$bandwidth
            neighbors <- nbh$neighbors
        }
        else{ # nearest neighbors within each stratum
            nbh.list <- lapply(split(Z,Sfactor),neighborhood,bandwidth=bandwidth)
            bandwidth <- sapply(nbh.list,function(nbh)nbh$bandwidth)
            tabS <- c(0,cumsum(tabulate(Sfactor))[-NS])
            neighbors <- unlist(lapply(1:NS,function(l){ ## incrementing the neighbors by
                nbh.list[[l]]$neighbors+tabS[l]}),use.names=FALSE) ## the size of the previous strata
        }
        response <- response[neighbors,,drop=FALSE]
    }
  
  # }}}
  # {{{ delay (left truncation)
  delayed <- attr(event.history,"entry.type")=="leftTruncated"
  ## && !(attr(event.history,"entry.type")=="")
  if (!delayed) { ## either NULL or ""
      entrytime <- NULL
  }  else {
      entrytime <- response[,"entry"]
      if(!(all(entrytime>=0)))
          stop(paste("Not all entry times in dataset are greater or equal to zero."))
  }
  # }}}

    # {{{  bound on the number of unique time points over all strata  
    switch(cotype,
           { # type=1
               size.strata <- NROW(response)
               NU <- 1
               if (cens.type!="intervalCensored")
                   N <- length(unique(response[,"time"]))
               else
                   N <- length(unique(response[,"L"]))
               ## if (delayed) N <- N + length(entrytime)
               if (delayed) N <- length(unique(c(entrytime,response[,"time"])))
           },
           {   # type=2
               size.strata <- tabulate(Sfactor)
               N <- NROW(response)
               NU <- length(size.strata)
               if (delayed) N <- 2*N
           },
           {   # type=3
               size.strata <- nbh$size.nbh
               N <- sum(size.strata)
               NU <- nbh$nu
               if (delayed) N <- 2*N
           },
           {   # type=4
               size.strata <- unlist(lapply(nbh.list,function(nbh)nbh$size.nbh),use.names=FALSE)
               N <- sum(size.strata)
               if (delayed) N <- 2*N
               n.unique.strata <- unlist(lapply(nbh.list,function(nbh)nbh$nu),use.names=FALSE)
               NU <- sum(n.unique.strata)
           })
  
  # }}}

  # {{{  characterizing the covariate space
  
    continuous.predictors <- colnames(NN)
    discrete.predictors <- colnames(strata)
    X <- switch(cotype,
                {#type=1
    NULL},
                { #type=2
    X <- data.frame(unique(strata[sorted,,drop=FALSE]))
    ## colnames(X) <- paste("strata",names(strata),sep=".")
    # colnames(X) <- names(strata)
    rownames(X) <- 1:NROW(X)
    X
},
                { #type=3
    X <- unlist(lapply(nbh.list,function(x)x$values),use.names=FALSE)
    X <- data.frame(X)
    ## colnames(X) <- paste("NN",names(NN),sep=".")
    colnames(X) <- colnames(NN)
    rownames(X) <- 1:NROW(X)                
    X
},
                { #type=4
    D <- data.frame(unique(strata[sorted,,drop=FALSE]))
    ## colnames(D) <- paste("strata",names(strata),sep=".")
    D <- data.frame(D[rep(1:NS,n.unique.strata),,drop=FALSE])
    C <- data.frame(unlist(lapply(nbh.list,function(x)x$values),use.names=FALSE))
    X <- cbind(D,C)
    ## colnames(X) <- c(paste("strata",names(strata),sep="."),paste("NN",names(NN),sep="."))
    colnames(X) <- c(colnames(strata),colnames(NN))
    rownames(X) <- 1:NROW(X)                
    X
},
                { #type=5
    X=data.frame(pseudo="pseudo")
    rownames(X) <- 1:NROW(X)                
    X
})
    if (x==TRUE)
        model.matrix <- switch(cotype,{NULL},strata,NN,cbind(strata,NN))[event.time.order,,drop=FALSE]
    else
        model.matrix <- NULL
    event.history <- event.history[event.time.order,,drop=FALSE]
    # }}}
    # {{{ caseweights
    if (missing(caseweights)) {
        weighted <- 0
        caseweights <- NULL
    }
    else {
        weighted <- 1
        if(length(caseweights)!=NROW(response))
            stop(paste("The length of caseweights is: ",
                       length(caseweights),
                       "\nthis is not the same as the number of subjects\nwith no missing values, which is ",
                       NROW(response),
                       sep=""))
        ## wrong to order by event.time.order when there are covariates
        ## caseweights <- caseweights[event.time.order]
        ## this fixes bug in versions < 1.5.7 
        caseweights <- caseweights[sorted]
    }
    # }}}
    # {{{  cluster correlated data need an adjusted variance formula
    clustered <- (length(covariates$cluster)>0)
    if (clustered)
        clustervar <- colnames(covariates$cluster)
    else
        clustervar <- NULL
    if (clustered){
        cluster <- covariates$cluster[sorted,,drop=TRUE]
        if (cotype==1){
            NC <- length(unique(cluster))
            cluster <- factor(cluster,labels=1:NC)
        }
        else{
            if (cotype==2){
                NC <- unlist(tapply(cluster,Sfactor,function(x){length(unique(x))}))
                cluster <- as.numeric(unlist(tapply(cluster,Sfactor,function(x){
                                                                factor(x,labels=1:length(unique(x)))})))
            }
        }
    }
    # }}}
    # {{{  find the appropriate C routine
    # with respect to model.type, cens.type, cotype and clustered
    # the following cases are not yet available
    ## if (length(attr(event.history,"entry.type"))>1) stop("Prodlim: Estimation for left-truncated data not yet implemented.")
    if (delayed & weighted>0) stop("Prodlim: Estimation for left-truncated data with caseweights not implemented.")
    if (reverse && cens.type!="rightCensored") stop("Prodlim: Estimation of the censoring distribution works only for right censored data.")
    if (delayed && clustered) stop("Prodlim: Estimation with delayed entry and cluster-correlated observations not yet implemented.")
    if (reverse && clustered) stop("Prodlim: Estimation of censoring distribution with cluster-correlated observations not yet handled.")
    if (cens.type=="intervalCensored" && model.type>=2) stop("Prodlim: Interval censored observations only handled for two-state models")
    ##   if (cens.type=="intervalCensored" && model.type>2) stop("Interval censored observations only handled for two-state and competing risks models")
    if (clustered && model.type>1) stop("Prodlim: Cluster-correlated observations only handled for two-state models")
    if (clustered && cotype %in% c(3,4)) stop("Prodlim: Cluster-correlated observations not yet handled in presence of continuous covariates") #cluster <- cluster[neighbors]
    if (cotype>1 && cens.type=="intervalCensored") stop("Prodlim: Interval censored data and covariate strata not yet handled.")
    if (model.type==1){
        # }}}
        # {{{  two state model
        if (clustered){
            ## right censored clustered
            fit <- .C("prodlim",as.double(response[,"time"]),as.double(response[,"status"]),integer(0),as.double(entrytime),as.double(caseweights),as.integer(cluster),as.integer(N),integer(0),as.integer(NC),as.integer(NU),as.integer(size.strata),time=double(N),nrisk=double(2*N),nevent=double(2*N),ncens=double(2*N),surv=double(N),cuminc=double(0),hazard=double(N),var.hazard=double(N+N),extra.double=double(4 * max(NC)),max.nc=as.integer(max(NC)),ntimes=integer(1),ntimes.strata=integer(NU),first.strata=integer(NU),reverse=integer(0),model=as.integer(0),independent=as.integer(0),delayed=as.integer(delayed),weighted=as.integer(weighted),PACKAGE="prodlim")
            NT <- fit$ntimes
            Cout <- list("time"=fit$time[1:NT],"n.risk"=matrix(fit$nrisk,ncol=2,byrow=FALSE,dimnames=list(NULL,c("n.risk","cluster.n.risk")))[1:NT,],"n.event"=matrix(fit$nevent,ncol=2,byrow=FALSE,dimnames=list(NULL,c("n.event","cluster.n.event")))[1:NT,],"n.lost"=matrix(fit$ncens,ncol=2,byrow=FALSE,dimnames=list(NULL,c("n.lost","cluster.n.lost")))[1:NT,],"surv"=fit$surv[1:NT],"se.surv"=fit$surv[1:NT]*sqrt(pmax(0,fit$var.hazard[N+(1:NT)])),"naive.se.surv"=fit$surv[1:NT]*sqrt(pmax(0,fit$var.hazard[1:NT])),"hazard"=fit$hazard[1:NT],"first.strata"=fit$first.strata,"size.strata"=fit$ntimes.strata,"model"="survival")
            Cout$maxtime <- max(Cout$time)
        }
        else{
            if (cens.type=="intervalCensored"){
                if (length(method)>1) method <- method[1]
                if (length(grep("impute",method))>0){
                    naiiveMethod <- strsplit(method,"impute.")[[1]][[2]]
                    if (naiiveMethod=="midpoint"){
                        naiveResponse <- data.frame(unclass(response))
                        naiveResponse$imputedTime <- (naiveResponse$L+naiveResponse$R)/2
                        naiveResponse[naiveResponse[,"status"]==0,"imputedTime"] <- naiveResponse[naiveResponse[,"status"]==0,"L"]
                        Cout <- prodlim(Hist(imputedTime,status!=0)~1,data=naiveResponse)
                        return(Cout)
                    }
                }
                else{
                    Cout <- prodlimIcensSurv(response,
                                             grid,
                                             tol=tol,
                                             maxiter=maxiter,
                                             ml=ifelse(method=="one.step",FALSE,TRUE),
                                             exact=exact)
                }
            }
            else{
                ## right censored not clustered
                fit <- .C("prodlim",as.double(response[,"time"]),as.double(response[,"status"]),integer(0),as.double(entrytime),as.double(caseweights),integer(0),as.integer(N),integer(0),integer(0),as.integer(NU),as.integer(size.strata),time=double(N),nrisk=double(N),nevent=double(N),ncens=double(N),surv=double(N),double(0),hazard = double(N),var.hazard=double(N),extra.double=double(0),max.nc=integer(0),ntimes=integer(1),ntimes.strata=integer(NU),first.strata=integer(NU),as.integer(reverse),model=as.integer(0),independent=as.integer(1),delayed=as.integer(delayed),weighted=as.integer(weighted),PACKAGE="prodlim")
                NT <- fit$ntimes
                Cout <- list("time"=fit$time[1:NT],
                             "n.risk"=fit$nrisk[1:NT],
                             "n.event"=fit$nevent[1:NT],
                             "n.lost"=fit$ncens[1:NT],
                             "surv"=fit$surv[1:NT],
                             "se.surv"=fit$surv[1:NT]*sqrt(pmax(0,fit$var.hazard[1:NT])),
                             "hazard"=fit$hazard[1:NT],
                             "first.strata"=fit$first.strata,
                             "size.strata"=fit$ntimes.strata,
                             "model"="survival")
                Cout$maxtime <- max(Cout$time)
            }
        }
    }
    else{
        # }}}
        # {{{ competing.risks model
        if (model.type==2){
            states <- attr(response,"states")
            E <- response[,"event"]-1 # for the c routine
            D <- response[,"status"]
            NS <- length(unique(E[D!=0])) # number of different causes
            fit <- .C("prodlim",
                      as.double(response[,"time"]),
                      as.double(D),
                      as.integer(E),
                      as.double(entrytime),
                      as.double(caseweights),
                      integer(0),
                      as.integer(N),
                      as.integer(NS),
                      integer(0),
                      as.integer(NU),
                      as.integer(size.strata),
                      time=double(N),
                      nrisk=double(N),
                      nevent=double(N * NS),
                      ncens=double(N),
                      surv=double(N),
                      cuminc=double(N * NS),
                      cause.hazard = double(N * NS),
                      var.hazard=double(N * NS),
                      extra.double=double(4 * NS),
                      max.nc=integer(0),
                      ntimes=integer(1),
                      ntimes.strata=integer(NU),
                      first.strata=integer(NU),
                      reverse=integer(0),
                      model=as.integer(1),
                      independent=as.integer(1),
                      delayed=as.integer(delayed),
                      weighted=as.integer(weighted),
                      PACKAGE="prodlim")
            NT <- fit$ntimes
            # changed Tue Sep 30 12:51:58 CEST 2008
            # its easier to work with a list than with a matrix
            #      gatherC <- function(x,dimR=fit$ntimes,dimC=NS,names=states){
            #        matrix(x[1:(dimR*dimC)],ncol=dimC,byrow=TRUE,dimnames=list(rep("",dimR),names))
            #      }
            gatherC <- function(x,dimR=fit$ntimes,dimC=NS,names=states){
                out <- split(x[1:(dimR*dimC)],rep(1:NS,dimR))
                names(out) <- names
                out
            }
            Cout <- list("time"=fit$time[1:NT],
                         "n.risk"=fit$nrisk[1:NT],
                         "n.event"=gatherC(fit$nevent),
                         "n.lost"=fit$ncens[1:NT],
                         "cuminc"=gatherC(fit$cuminc),
                         "var.cuminc"=gatherC(fit$var.hazard),
                         "se.cuminc"=gatherC(sqrt(pmax(0,fit$var.hazard))),
                         "surv"=fit$surv[1:NT],
                         "cause.hazard"=gatherC(fit$cause.hazard),
                         "first.strata"=fit$first.strata,
                         "size.strata"=fit$ntimes.strata,
                         "model"="competing.risks")
            Cout$maxtime <- max(Cout$time)
        }
        else {
            # multi.state model
            # --------------------------------------------------------------------    
            Cout <- prodlimMulti(response,size.strata,N,NU)
            Cout$maxtime <- max(Cout$time)
        }
    }
    if (conf.int==TRUE) conf.int <- 0.95
    # }}}
    # {{{  confidence intervals
    if (is.numeric(conf.int) && cens.type!="intervalCensored"){
        if (model.type==1){
            if (!(is.null(Cout$se.surv))){
                ## pointwise confidence intervals for survival probability
                zval <- qnorm(1- (1-conf.int)/2, 0,1)
                lower <- pmax(Cout$surv - zval * Cout$se.surv,0)
                lower[Cout$se.surv==0] <- 0
                upper <- pmin(Cout$surv + zval * Cout$se.surv,1)
                upper[Cout$se.surv==0] <- 1
                Cout <- c(Cout,list(lower=lower,upper=upper))
            }
        }
        else{
            if (is.numeric(conf.int)){
                if (!(0<conf.int && conf.int<1)) conf.int <- 0.95
                ## pointwise confidence intervals for cumulative incidence probabilities
                # variance for cuminc (Korn & Dorey (1992), Stat in Med, Vol 11, page 815)
                zval <- qnorm(1- (1-conf.int)/2, 0,1)
                lower <- lapply(1:NS, function(state){
                                    pmax(Cout$cuminc[[state]] - zval * Cout$se.cuminc[[state]],0)})
                upper <- lapply(1:NS, function(state){
                                    pmin(Cout$cuminc[[state]] + zval * Cout$se.cuminc[[state]],1)})
                names(lower) <- states
                names(upper) <- states
                Cout <- c(Cout,list(lower=lower,upper=upper))
            }
        }
    }
    # }}}
    # {{{ return object of class "prodlim"
    out <- list("call"=call,
                "formula"=formula,
                "model.response"=event.history,
                "originalDataOrder"=order(event.time.order),
                "X"=X,
                "model.matrix"=model.matrix,
                "discrete.predictors"=discrete.predictors,
                "continuous.predictors"=continuous.predictors,
                "xlevels"=xlevels,
                "clustervar"=clustervar,
                "covariate.type"=cotype,
                "cens.type"=cens.type,
                "conf.int"=conf.int,
                "reverse"=reverse,
                "type"=type,
                "na.action"=attr(EHF,"na.action"))
    if (cotype %in% c(3,4)) out <- c(out,list("bandwidth"=bandwidth))
    out <- c(Cout,out)
    class(out) <-  "prodlim"
    return(out)
    # }}}
}
