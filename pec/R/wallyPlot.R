# {{{ Wally plot function
##' Wally plots to assess calibration of a risk or survival prediction
##'
##' @title Wally plots to assess calibration of a risk or survival prediction
##' @param object Probabilistic survival predictions or probabilistic event risk predictions
##' evaluated at \code{time} for the subjects in \code{data}. Either
##' given in form of a numeric vector of probabilistic predictions or
##' as an object which in the survival setting has a
##' \code{predictSurvProb} method and in the competing risks setting
##' has a \code{predictEventProb} method.
##' @param time Time interest for evaluating calibration of the
##' predictions.
##' @param formula A survival or event history formula. The left hand
##' side is used to compute the expected event status. If
##' \code{formula} is \code{missing}, try to extract a formula from
##' the first element in object.
##' @param data A data frame in which to validate the prediction
##' models and to fit the censoring model. If \code{data} is missing,
##' try to extract a data set from the first element in object.
##' @param cause For competing risks settings the cause of interest.
##' @param q The number of quantiles. Defaults to 10.
##' @param ylim Limits of y-axis. If missing the function tries to
##' find appropriate limits based on the simulated and real data.
##' @param hanging If \code{TRUE}, hang bars corresponding to observed
##' frequencies at the value of the corresponding prediction.
##' @param seed A seed value to make results reproducible.
##' @param mar Plot margins passed to par.
##' @param colbox Color of the box which identifies the real data
##' calibration plot.
##' @param type For survival models only: show either "risk" or
##' "survival".
##' @param pseudo Logical. Determines the method for estimating expected event frequencies. See \code{calPlot}. Default is \code{FALSE}.
##' @param verbose If \code{TRUE} warn about missing formula and data.
##' @param ... Further arguments passed to \code{calPlot}.
##' @return List of simulated and real data.
##' @seealso calPlot
##' @examples
##'
##' # Survival setting
##' library(prodlim)
##' library(survival)
##' set.seed(180)
##' d = SimSurv(180)
##' f = coxph(Surv(time,status)~X1+X2,data=d)
##' \dontrun{
##' wallyPlot(f,
##'           time=4,
##'           q=10,
##'           type="risk",
##'           data=d,
##'           formula=Surv(time,status)~1)
##'  wallyPlot(f,
##'           time=4,
##'           q=10,
##'           hanging=TRUE,
##'           type="survival",
##'           data=d,
##'           formula=Surv(time,status)~1)
##' }
##' 
##' # Competing risks setting
##' library(prodlim)
##' library(survival)
##' library(riskRegression)
##' set.seed(180)
##' d2 = SimCompRisk(180)
##' f2 = CSC(Hist(time,event)~X1+X2,data=d2)
##' \dontrun{
##' wallyPlot(f2,
##'           time=5,
##'           q=3,
##'           hanging=TRUE,
##'           data=d2,
##'           formula=Hist(time,event)~1)
##'           
##' }
##' 
##' @export 
##' @author Paul F. Blanche <paul.blanche@@univ-ubs.fr> and Thomas A. Gerds <tag@@biostat.ku.dk>
wallyPlot <- function(object,
                      time,
                      formula,
                      data,
                      cause=1,
                      q=10,
                      ylim,
                      hanging=FALSE,
                      seed=NULL,
                      mar=c(4.1,4.1, 2, 2),
                      colbox="red",
                      type="risk",
                      pseudo=FALSE,
                      ## identify="select from list",
                      verbose=TRUE,...){
    identify="select from list"
    type <- match.arg(type,c("risk","survival"))
    # {{{ data & formula
    if (missing(data)){
        trydata <- try(data <- eval(object$call$data),silent=TRUE)
        if (("try-error" %in% class(trydata))|| match("data.frame",class(data),nomatch=0)==0)
            stop("Argument data is missing.")
        else
            if (verbose)
                warning("Argument data is missing. I use the data from the call to the first model instead.")
    }
    if (missing(formula)){
        tryformula <- try(as.character(object$call$formula),silent=TRUE)
        if (("try-error" %in% class(tryformula))||length(grep("~",as.character(object$call$formula)))==0){
            stop(paste("Argument formula is missing and first model has no usable formula."))
        } else{
              ftry <- try(formula <- eval(object$call$formula),silent=TRUE)
              if ((class(ftry)=="try-error") || match("formula",class(formula),nomatch=0)==0)
                  stop("Argument formula is missing and first model has no usable formula.")
              else if (verbose)
                  warning("Formula missing. Using formula from first model")
              ## remove covariates
              formula <- update.formula(formula,".~1")
          }
    }
    m <- model.frame(formula,data,na.action=na.fail)
    response <- model.response(m)
    if (match("Surv",class(response),nomatch=FALSE))
        model.type <- "survival"
    else
        model.type <- attr(response,"model")
    if (is.null(model.type) & length(unique(response))==2)
        model.type <- "binary"
    if (!(model.type=="binary")){
        neworder <- order(response[,"time"],-response[,"status"])
        response <- response[neworder,,drop=FALSE]
        Y <- response[,"time"]
        status <- response[,"status"]
        if (model.type=="competing.risks"){
            event <- as.numeric(response[,"event"])
            event[status==0] <- 0
            status <- event
            rm(event)
        }
        data <- data[neworder,]
        if (missing(time))
            time <- median(Y)
        else
            if (length(time)>1)
                stop("Please specify only one time point.")
    }
    predictHandlerFun <- switch(model.type,
                                "binary"="predictStatusProb",
                                "competing.risks"="predictEventProb",
                                "survival"="predictSurvProb")
    # }}}
    # {{{ arguments for calPlot
    calPlot.defaultArgs <- list(formula=formula,
                                cause=cause,
                                time=time,
                                q=q,
                                type=type,
                                na.action=na.fail,
                                col=c("grey90","grey30"),
                                ylim=c(0,1),
                                percent=TRUE,
                                pseudo=pseudo,
                                showPseudo=FALSE,
                                bars=TRUE,
                                names=FALSE,
                                hanging=hanging,
                                showFrequencies=FALSE,
                                legend=FALSE,
                                axis2.mgp=c(4,1,0),
                                method="quantile",
                                q=10,
                                ylab="",
                                xlab="")
    superuser.defaults <- list(hide=TRUE,choice=6,zoom=FALSE)
    smartA <- prodlim::SmartControl(call= list(...),
                                    keys=c("calPlot","superuser"),
                                    ignore=NULL,
                                    ignore.case=TRUE,
                                    defaults=list("calPlot"=calPlot.defaultArgs,"superuser"=superuser.defaults),
                                    verbose=TRUE)
    # }}}
    # {{{ define risk groups
    if (class(object)[1] %in% c("numeric","double"))
        object <- matrix(object,ncol=1)
    pred <- switch(model.type,
                   "competing.risks"={
                       p <- as.vector(do.call(predictHandlerFun,list(object,newdata=data,times=time,cause=cause)))
                       if (class(object)[[1]]%in% c("matrix","numeric")) p <- p[neworder]
                       p
                   },
                   "survival"={
                       p <- as.vector(do.call(predictHandlerFun,list(object,newdata=data,times=time)))
                       if (class(object)[[1]]%in% c("matrix","numeric")) p <- p[neworder]
                       p ## the prediction has to be a survival probability!
                   },
                   "binary"={
                       p <- do.call(predictHandlerFun,list(object,newdata=data))
                       if (class(object)[[1]]%in% c("matrix","numeric")) p <- p[neworder]
                       p
                   })
    if (any(is.na(pred))) stop("Missing values in prediction. Maybe time interest needs to be set earlier?")
    quant <- quantile(pred,seq(0,1,1/q))
    Pred.cut <- cut(pred,breaks =quant ,labels=1:(length(quant)-1),include.lowest=TRUE)
    ## Test if time is ok
    if (any((tooshortFollowup <- tapply(Y,Pred.cut,max))<time))
        stop(paste0("The following risk groups have too short followup. Nth quantile of predicted risks: ",paste(names(tooshortFollowup),collapse=", "),"\nThe minimum of the maximal followup in the risk groups is: ",
                    signif(min(tooshortFollowup),2)))
    if ("data.table" %in% class(data))
        MyData <- cbind(data[,all.vars(update(formula,".~1")),drop=FALSE,with=FALSE],Pred.cut=Pred.cut,pred=pred)
    else
        MyData <- cbind(data[,all.vars(update(formula,".~1")),drop=FALSE],Pred.cut=Pred.cut,pred=pred)
    ## MyData <- data.frame(time=Y,status=status,Pred.cut=Pred.cut,pred=pred)
    form.pcut <- update(formula,paste(".~Pred.cut"))
    if (model.type=="competing.risks"){
        # estimate the cumulative incidence of the competing event at time t in each subgroup
        crFit <- prodlim::prodlim(form.pcut,data=MyData)
        # create the vector of prediction for cause 2, based on the sub-group specific cumulative incidence
        pred2 <- unlist(predict(crFit,
                                cause=2,
                                times=time,
                                newdata=data.frame(Pred.cut=Pred.cut),
                                type="cuminc"))
    }
    # }}}
    # {{{ functions to draw event times
    GenTimesFromPred <- function(pi,t){
        # compute the parameters for the exponential distributions
        lambdai <- -(1/t)*log(1-pmin(pi,1))
        Ti <- rexp(n=length(pi),rate=lambdai)
        Ti
    }
    GenCensFromKMfit <- function(times,status,groups=NULL){
        if (all(status!=0)){
            rep(Inf,length(times))
        }else{
             if (is.null(groups)){
                 # uniform numbers
                 U <- runif(length(times),0,1)
                 # fit marginal KM for censoring
                 fitcens <- prodlim::prodlim(Hist(time,status)~1,data=data.frame(time=times,status=status),reverse=TRUE)
                 # max time in data
                 tau <- max(fitcens$time)
                 c(fitcens$time,tau)[prodlim::sindex(jump.times=c(0,1-fitcens$surv),eval.times=U)]
             } else{
                   # fit stratified KM for censoring
                   fitcens <- prodlim::prodlim(Hist(time,status)~groups,data=data.frame(time=times,
                                                                            status=status,
                                                                            groups=groups),reverse=TRUE)
                   n.groups <- table(groups)
                   C <- unlist(lapply(1:length(n.groups),function(g){
                                                  U.g <- runif(n.groups[g],0,1)
                                                  start <- fitcens$first.strata[g]
                                                  size <- fitcens$size.strata[g]
                                                  strata.g <- start:(start+size)
                                                  g.time <- fitcens$time[strata.g]
                                                  tau.g <- max(g.time,na.rm=TRUE)
                                                  g.surv <- fitcens$surv[strata.g]
                                                  isna <- is.na(g.time)
                                                  g.surv[isna] <- 0
                                                  g.time[isna] <- tau.g
                                                  cgroup <- c(g.time,tau.g)[prodlim::sindex(jump.times=c(0,1-g.surv),eval.times=U.g)]
                                                  cgroup
                                              }))
                   C[order(order(groups))]
               }
         }
    }
    GetResponseSurv <- function(Y,status,pred,Pred.cut,time){
        C <- GenCensFromKMfit(times=Y,status=status,groups=Pred.cut)
        ## Note: predictions pred are survival probabilities
        ##       we simulate event times based on 1-pred
        T <- GenTimesFromPred(pi=1-pred,t=time)
        Tcens <- pmin(T,C)
        delta <- as.numeric(T<=C)
        data.frame(time=Tcens,status=delta,pred=pred,Pred.cut=Pred.cut)
    }
    GetResponseCR <- function(Y,status,pred,pred2,Pred.cut,time){
        # generate censoring conditionally on groups
        C <- GenCensFromKMfit(times=Y,status=status,groups=Pred.cut)
        T <- GenTimesFromPred(pi=(pred+pred2),t=time)
        Tcens <- pmin(T,C)
        lambda1 <- -1/( time*(1+pred2/pred) )*log(1-pmin(pred2+pred,0.999999))
        lambda2 <- lambda1*(pred2/pred)
        eta <- rbinom(length(T),1,prob=lambda2/(lambda1+lambda2))+rep(1,length(T))               
        delta <- as.numeric(T<C)*eta
        data.frame(time=Tcens,status=delta,pred=pred,Pred.cut=Pred.cut)
    }
    # }}}
    # set the seed for making the Wally plot repoducible
    set.seed(seed)
    # {{{ loop creates the simulated data sets
    if(model.type!="competing.risks"){
        DataList <- lapply(1:8,function(i){
                               GetResponseSurv(Y=Y,
                                               status=status,
                                               pred=pred, 
                                               Pred.cut=Pred.cut,
                                               time=time)
                           })
    }else{
         DataList <- lapply(1:8,function(i){
                                GetResponseCR(Y=Y,
                                              status=status,
                                              pred=pred,
                                              pred2=pred2,
                                              Pred.cut=Pred.cut,
                                              time=time)
                            })
     }
    ## add real data at the last position
    DataList <- c(DataList,list(MyData))
    # }}}  
    # {{{ sample the order of the data in the list
    pos <- 1:9
    pos <- sample(pos)        
    figpos <- order(pos)[9] # where is the actual plot (with the true data)
    # }}}
    # {{{ prepare graphics output
    oldpar <- par(no.readonly = TRUE)
    par(mar = mar)
    par(mfrow=c(3,3),oma=c(3,0,0,0))
    # }}}
    # {{{ loop to create the 9 plots
    TabList <- vector("list",9)
    printleg <- c(rep(FALSE,4),TRUE,rep(FALSE,4))
    for(i in pos){
        smartA$calPlot$data <- DataList[[i]]
        if (i==9) smartA$calPlot$formula <- formula
        else smartA$calPlot$formula <- formula("Hist(time,status)~1")
        ## if (model.type=="survival")
        ## smartA$calPlot$object <- 1- DataList[[i]]$pred
        ## else            
        smartA$calPlot$object <- DataList[[i]]$pred
        smartA$calPlot$plot <- FALSE
        TabList[[i]] <- do.call("calPlot",smartA$calPlot)
        if (any(is.na(TabList[[i]]$plotFrame[[1]]$Obs)))
            stop("Missing values in expected frequencies. Maybe too many quantiles relative to the number of observations? Or time interest set too late?")
    }
    if (hanging){
        rangeY <- range(sapply(TabList,function(x){
                                   c(min(x$plotFrame[[1]]$Pred-x$plotFrame[[1]]$Obs),max(x$plotFrame[[1]]$Pred))}))
        minY <- min(0,seq(-1,1,0.05)[prodlim::sindex(eval.times=rangeY[1],jump.times=seq(-1,1,0.05))])
    } else{
          if (model.type=="survival" && type=="risk") 
              rangeY <- range(sapply(TabList,function(x){range(c(1-x$plotFrame[[1]]))}))
          else
              rangeY <- range(sapply(TabList,function(x){range(c(x$plotFrame[[1]]))}))
          minY <- 0
      }
    maxY <- c(seq(0,1,0.05),1)[1+prodlim::sindex(eval.times=rangeY[2],jump.times=seq(0,1,0.05))]
    if (!missing(ylim)){
        if (minY<0) minY <- min(ylim,minY)
        maxY <- ylim[2]
    }
    for (j in 1:9){
        i = pos[j]
        px <- TabList[[i]]
        px$control$barplot$ylim <- c(minY,maxY)
        ## need to round to avoid strange results like in
        ## paste(100*seq(-0.15,0.7,0.7/4),"%")
        ## px$control$axis2$at <- sort(c(0,round(seq(minY,maxY,(maxY-minY)/4),3)))
        px$control$axis2$at <- sort(round(seq(minY,maxY,(maxY-minY)/4),3))
        if (j==8) {
            px$legend=TRUE
            if (model.type=="survival" && type=="survival"){
                px$control$barplot$legend.text=c("Predicted survival","Observed frequency")
                px$control$legend$legend <- c("Predicted survival","Observed frequency")
            }
            else{
                px$control$barplot$legend.text=c("Predicted risk","Observed frequency")
                px$control$legend$legend <- c("Predicted risk","Observed frequency")
            }
            px$control$legend$xpd <- NA
            px$control$legend$x <- "bottom"
            ## px$control$legend$x <- 0
            px$control$legend$ncol=2
            px$control$legend$inset=c(0,-0.4)
            ## px$control$legend$y <- minY - (par()$mai[1]/par()$mar[1])/2
            px$control$legend$cex <- 1.3
        }
        ## add the plot to the grid
        plot(px)
        if (identify!="click"){
            upleft <- par("usr")[c(1,4)]
            points(x=upleft[1],y=upleft[2],xpd=NA,pch=19,cex=5,col="orange")
            text(j,x=upleft[1],y=upleft[2],xpd=NA,col="black",cex=1.2)
        }
    }
    ## mtext(side=1,line=1,"Can you find wally?",cex=1.5*par()$cex)
    # }}}
    # {{{ ask to show the actual plot
    if (is.null(smartA$superuser$choice)){
        par(oldpar)
        invisible(TabList[order(pos)])
    }else{
         if (smartA$superuser$hide!=FALSE){
             if (identify=="click"){
                 stop("Does not work yet")
                 ## par(mfg = c(1,1))
                 ## ans <- identify(x=c(0,1),y=NULL)
                 ## usr <- par("usr")
                 ## xx <- 1
                 ## number <- par$
                 ## identify(x
             }else{
                  xx <-  select.list(1:9,
                                     multiple=FALSE,
                                     title="\nWhere is Wally? Can you find the plot which is based on the real data?\nSelect an orange number: ")
              }
         }else xx <- smartA$superuser$choice
         par(mfg = c(xx%/%3.1 + 1, xx - (xx%/%3.1) * 3))
         box(col = "green", lwd = 3)
         # }}}
         # {{{ show the actual plot by adding a red box
         if (xx==figpos) {
             mtext("Correct: real data!",col="green",cex=1.2,line=-3,xpd=NA)
         }else{
              mtext("Not correct!",col="red",cex=1.2,line=-3,xpd=NA)
              par(mfg = c(figpos%/%3.1 + 1, figpos - (figpos%/%3.1) * 3))
              mtext("Real data",col="red",cex=1.2,line=-3,side=3,xpd=NA)
              box(col = colbox, lwd = 5)
          }
         # }}}
         # {{{ ask to press a key to zoom in on the actual plot
         ## readline("Hit <Enter> to better show the original plot. ")
         if (smartA$superuser$zoom==FALSE && smartA$superuser$hide!=FALSE){
             zoom <- select.list(c("yes","no"),title="Zoom in on real data calibration plot? ")
         }else zoom <- "no"
         # }}}
         # {{{ show the actual plot
         if((smartA$superuser$zoom==TRUE) || (zoom=="yes")){
             par(mfrow=c(1,1),oma=c(2,2,2,2),mar=c(4.1,4.1, 4.1, 4.1))
             px <- TabList[[9]]
             px$control$barplot$ylim <- c(minY,maxY)
             ## need to round to avoid strange results like in
             ## paste(100*seq(-0.15,0.7,0.7/4),"%")
             ## px$control$axis2$at <- sort(c(0,round(seq(minY,maxY,(maxY-minY)/4),3)))
             px$control$axis2$at <- sort(round(seq(minY,maxY,(maxY-minY)/4),3))
             px$legend=TRUE
             if (model.type=="survival" && type=="survival"){
                 px$control$barplot$legend.text=c("Predicted survival","Observed frequency")
                 px$control$legend$legend <- c("Predicted survival","Observed frequency")
                 px$control$barplot$xlab <- "Survival groups"
             }else{
                  px$control$barplot$legend.text=c("Predicted risk","Observed frequency")
                  px$control$legend$legend <- c("Predicted risk","Observed frequency")
                  px$control$barplot$xlab <- "Risk groups"
              }
             ## print(px$control$barplot$legend.text)
             px$control$legend$xpd <- NA
             px$control$legend$x <- "top"
             px$control$legend$ncol <- 2
             px$control$legend$inset <- c(0,-0.2)
             ## px$control$legend$cex <- 1.3*par()$cex
             ## px$control$names$cex <- par()$cex
             ## px$control$frequencies$cex <- par()$cex
             px$showFrequencies <- TRUE
             qq <- attr(px$plotFrames[[1]],"quantiles")
             px$names <- paste0(sprintf("%1.1f",100*qq[-length(qq)])," - ",
                                sprintf("%1.1f",100*qq[-1]))
             plot(px)
         }
         # }}}
         par(oldpar)
         invisible(TabList[order(pos)])
         ## invisible(DataList)
     }
}
# }}}
