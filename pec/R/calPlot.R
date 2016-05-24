#' Calibration plots for right censored data
#' 
#' Calibration plots for risk prediction models in right censored survival and
#' competing risks data
#' 
#' For method "nne" the optimal bandwidth with respect to is obtained with the
#' function \code{\link{dpik}} from the package \code{KernSmooth} for a box
#' kernel function.
#' 
#' @param object A named list of prediction models, where allowed
#' entries are (1) R-objects for which a \link{predictSurvProb} method
#' exists (see details), (2) a \code{call} that evaluates to such an
#' R-object (see examples), (3) a matrix with predicted probabilities
#' having as many rows as \code{data} and as many columns as
#' \code{times}. For cross-validation all objects in this list must
#' include their \code{call}.
#' @param time The evaluation time point at predicted event
#' probabilities are plotted against pseudo-observed event status.
#' @param formula A survival or event history formula. The left hand
#' side is used to compute the expected event status. If
#' \code{formula} is \code{missing}, try to extract a formula from the
#' first element in object.
#' @param data A data frame in which to validate the prediction models
#' and to fit the censoring model. If \code{data} is missing, try to
#' extract a data set from the first element in object.
#' @param splitMethod Defines the internal validation design:
#' 
#' \code{none/noPlan}: Assess the models in the give \code{data}, usually
#' either in the same data where they are fitted, or in independent test data.
#' 
#' \code{BootCv}: Bootstrap cross validation. The prediction models
#' are trained on \code{B} bootstrap samples, that are either drawn
#' with replacement of the same size as the original data or without
#' replacement from \code{data} of the size \code{M}.  The models are
#' assessed in the observations that are NOT in the bootstrap sample.
#' @param B The number of cross-validation steps.
#' @param M The size of the subsamples for cross-validation.
#' @param pseudo Logical. Determines the method for estimating expected event status:
#' 
#' \code{TRUE}: Use average pseudo-values.  \code{FALSE}: Use
#' the product-limit estimate, i.e., apply the Kaplan-Meier method for
#' right censored survival and the Aalen-Johansen method for right
#' censored competing risks data. 
#' @param type Either "risk" or "survival". 
#' @param showPseudo If \code{TRUE} the
#' pseudo-values are shown as dots on the plot (only when \code{pseudo=TRUE}).
#' @param pseudo.col Colour for pseudo-values.
#' @param pseudo.pch Dot type (see par) for pseudo-values.
#' @param method The method for estimating the calibration curve(s):
#' 
#' \code{"nne"}: The expected event status is obtained in the nearest
#' neighborhood around the predicted event probabilities.
#' 
#' \code{"quantile"}: The expected event status is obtained in groups
#' defined by quantiles of the predicted event probabilities.
#' @param round If \code{TRUE} predicted probabilities are rounded to
#' two digits before smoothing. This may have a considerable effect on
#' computing efficiency in large data sets.
#' @param bandwidth The bandwidth for \code{method="nne"}
#' @param q The number of quantiles for \code{method="quantile"} and \code{bars=TRUE}.
#' @param bars If \code{TRUE}, use barplots to show calibration.
#' @param hanging  Barplots only. If \code{TRUE}, hang bars corresponding to observed frequencies at the value of the corresponding prediction.
#' @param names Barplots only. Names argument passed to \code{names.arg} of \code{barplot}.
#' @param showFrequencies Barplots only. If \code{TRUE}, show frequencies above the bars.
#' @param jack.density Gray scale for pseudo-observations.
#' @param plot If \code{FALSE}, do not plot the results, just return a plottable object.
#' @param add If \code{TRUE} the line(s) are added to an existing
#' plot.
#' @param diag If \code{FALSE} no diagonal line is drawn.
#' @param legend If \code{FALSE} no legend is drawn.
#' @param axes If \code{FALSE} no axes are drawn.
#' @param xlim Limits of x-axis.
#' @param ylim Limits of y-axis.
#' @param xlab Label for y-axis.
#' @param ylab Label for x-axis.
#' @param col Vector with colors, one for each element of
#' object. Passed to \code{\link{lines}}.
#' @param lwd Vector with line widths, one for each element of
#' object. Passed to \code{\link{lines}}.
#' @param lty lwd Vector with line style, one for each element of
#' object.  Passed to \code{\link{lines}}.
#' @param pch Passed to \code{\link{points}}.
#' @param cause For competing risks models, the cause of failure or
#' event of interest
#' @param percent If TRUE axes labels are multiplied by 100 and thus
#' interpretable on a percent scale.
#' @param giveToModel List of with exactly one entry for each entry in
#' \code{object}. Each entry names parts of the value of the fitted
#' models that should be extracted and added to the value.
#' @param na.action Passed to \code{\link{model.frame}}
#' @param cores Number of cores for parallel computing.  Passed as
#' value of argument \code{mc.cores} to \code{\link{mclapply}}.
#' @param verbose if \code{TRUE} report details of the progress,
#' e.g. count the steps in cross-validation.
#' @param cex Default cex used for legend and labels.
#' @param ... Used to control the subroutines: plot, axis, lines, barplot,
#' legend. See \code{\link{SmartControl}}.
#' @return list with elements: time, pseudoFrame and bandwidth (NULL for method
#' quantile).
#' @keywords survival 
##' @examples
##' 
##' library(prodlim)
##' library(lava)
##' library(riskRegression)
##' library(survival)
##' # survival
##' dlearn <- SimSurv(40)
##' dval <- SimSurv(100)
##' f <- coxph(Surv(time,status)~X1+X2,data=dlearn)
##' cf=calPlot(f,time=3,data=dval)
##' print(cf)
##' plot(cf)
##' 
##' g <- coxph(Surv(time,status)~X2,data=dlearn)
##' cf2=calPlot(list("Cox regression X1+X2"=f,"Cox regression X2"=g),
##'     time=3,
##'     type="risk",
##'     data=dval)
##' print(cf2)
##' plot(cf2)
##' calPlot(f,time=3,data=dval,type="survival")
##' calPlot(f,time=3,data=dval,bars=TRUE,pseudo=FALSE)
##' calPlot(f,time=3,data=dval,bars=TRUE,type="risk",pseudo=FALSE)
##' 
##' calPlot(f,time=3,data=dval,bars=TRUE,hanging=TRUE)
##' calPlot(f,time=3,data=dval,bars=TRUE,type="risk",hanging=TRUE)
##' 
##' set.seed(13)
##' m <- crModel()
##' regression(m, from = "X1", to = "eventtime1") <- 1
##' regression(m, from = "X2", to = "eventtime1") <- 1
##' m <- addvar(m,c("X3","X4","X5"))
##' distribution(m, "X1") <- binomial.lvm()
##' distribution(m, "X4") <- binomial.lvm()
##' d1 <- sim(m,100)
##' d2 <- sim(m,100)
##' csc <- CSC(Hist(time,event)~X1+X2+X3+X4+X5,data=d1)
##' fgr <- FGR(Hist(time,event)~X1+X2+X3+X4+X5,data=d1,cause=1)
##' predict.crr <- cmprsk:::predict.crr
##' cf3=calPlot(list("Cause-specific Cox"=csc,"Fine-Gray"=fgr),
##'         time=5,
##'         legend.x=-0.3,
##'         legend.y=1.35,
##'         ylab="Observed event status",
##'         legend.legend=c("Cause-specific Cox regression","Fine-Gray regression"),
##'         legend.xpd=NA)
##' print(cf3)
##' plot(cf3)
##' 
##' b1 <- calPlot(list("Fine-Gray"=fgr),time=5,bars=TRUE,hanging=FALSE)
##' print(b1)
##' plot(b1)
##' 
##' calPlot(fgr,time=5,bars=TRUE,hanging=TRUE)
##' 
##' 
#' @author Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#' @export 
calPlot <- function(object,
                    time,
                    formula,
                    data,
                    splitMethod="none",
                    B=1,
                    M,
                    pseudo,
                    type,
                    showPseudo,
                    pseudo.col=NULL,
                    pseudo.pch=NULL,
                    method="nne",
                    round=TRUE,
                    bandwidth=NULL,
                    q=10,
                    bars=FALSE,
                    hanging=FALSE,
                    names="quantiles",
                    showFrequencies=FALSE,
                    jack.density=55,
                    plot=TRUE,
                    add=FALSE,
                    diag=!add,
                    legend=!add,
                    axes=!add,
                    xlim=c(0,1),
                    ylim=c(0,1),
                    xlab,
                    ylab,
                    col,
                    lwd,
                    lty,
                    pch,
                    cause=1,
                    percent=TRUE,
                    giveToModel=NULL,
                    na.action=na.fail,
                    cores=1,
                    verbose=FALSE,
                    cex=1,
                    ...){
    if (missing(pseudo)){
        if(method=="quantiles"||bars==TRUE) 
            pseudo <- FALSE
        else pseudo <- TRUE
    }
    if (missing(showPseudo))
        showPseudo <- ifelse(add||(pseudo!=FALSE),FALSE,TRUE)
    # {{{ find number of objects and lines
    cobj=class(object)[[1]]
    if (cobj!="list"){
        object <- list(object)
    }
    if (is.null(names(object))) names(object) <- paste0("Model.",1:length(object))
    if (bars){
        method="quantile"
        if (!(length(object)==1)) stop(paste0("Barplots work only for one prediction at a time. Provided are ",length(object), "predictions"))
    }
    if (is.null(names(object))){
        names(object) <- sapply(object,function(o)class(o)[1])
        names(object) <- make.names(names(object),unique=TRUE)
    }
    else{
        names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
    }
    NF <- length(object)
    # }}}
    # {{{ lines types
    if (missing(lwd)) lwd <- rep(3,NF)
    if (missing(col)) {
        if (bars)
            col <- c("grey90","grey30")
        else
            col <- 1:NF
    }
    if (missing(lty)) lty <- rep(1, NF)
    if (missing(pch)) pch <- rep(1, NF)
    if (length(lwd) < NF) lwd <- rep(lwd, NF)
    if (length(lty) < NF) lty <- rep(lty, NF)
    if (length(col) < NF) col <- rep(col, NF)
    if (length(pch) < NF) pch <- rep(pch, NF)
    # }}}
    # {{{ data & formula
    if (missing(data)){
        data <- eval(object[[1]]$call$data)
        if (match("data.frame",class(data),nomatch=0)==0)
            stop("Argument data is missing.")
        else
            if (verbose)
                warning("Argument data is missing. I use the data from the call to the first model instead.")
    }
    if (missing(formula)){
        if (length(grep("~",as.character(object[[1]]$call$formula)))==0){
            stop(paste("Argument formula is missing and first model has no usable formula:",as.character(object[[1]]$call$formula)))
        } else{
              ftry <- try(formula <- eval(object[[1]]$call$formula),silent=TRUE)
              if ((class(ftry)=="try-error") || match("formula",class(formula),nomatch=0)==0)
                  stop("Argument formula is missing and first model has no usable formula.")
              else if (verbose)
                  warning("Formula missing. Using formula from first model")
          }
    }
    
    m <- model.frame(formula,data,na.action=na.action)
    response <- model.response(m)
    if (match("Surv",class(response),nomatch=FALSE))
        model.type <- "survival"
    else
        model.type <- attr(response,"model")
    if (is.null(model.type) & length(unique(response))==2)
        stop("This function works only for survival and competing risks models.")
    ## model.type <- "binary"
    if (missing(type))
        type <- ifelse(model.type=="survival","survival","risk")
    if (missing(ylab))
        if (bars)
            ylab=""
        else
            ylab <- ifelse(type=="survival","Observed survival frequencies","Observed event frequencies")
    if (type=="survival" && !(model.type %in% c("survival","binary")))
        stop(paste0("Type survival works only in survival or binary outcome models. This is a ",model.type, " model"))
    if (!(model.type=="binary")){
        neworder <- order(response[,"time"],-response[,"status"])
        response <- response[neworder,,drop=FALSE]
        Y <- response[,"time"]
        ## status <- response[,"status"]
        data <- data[neworder,]
        # }}}
        # {{{ prediction timepoint 

        if (missing(time))
            time <- median(Y)
        else
            if (length(time)>1)
                stop("Please specify only one time point.")
    }

    # }}}
    # {{{ compute pseudo-values

    #  require(pseudo)
    #  jack=pseudosurv(time=Y,event=status,tmax=time)[[3]]
    predictHandlerFun <- switch(model.type,
                                "binary"="predictStatusProb",
                                "competing.risks"="predictEventProb",
                                "survival"="predictSurvProb")
    if (pseudo==FALSE && splitMethod!="none")
        stop(paste0("Split method ",splitMethod," is only implemented for : 'pseudo=TRUE'."))
    if (model.type=="binary")
        if (is.factor(response))
            jack <- as.numeric(response==levels(response)[2])
        else
            jack <- as.numeric(response)
    ## ==levels(response)[1])
    else{
        if (pseudo==TRUE){
            margForm <- update(formula,paste(".~1"))
            margFit <- prodlim::prodlim(margForm,data=data)
            jack <- prodlim::jackknife(margFit,cause=cause,times=time)
        }else{## prodlim in strata defined by predictions
             jack <- NULL
         }
    }
    # }}}
    # {{{ smartControl
    axis1.DefaultArgs <- list(side=1,las=1,at=seq(0,ylim[2],ylim[2]/4))
    axis2.DefaultArgs <- list(side=2,las=2,at=seq(0,ylim[2],ylim[2]/4),mgp=c(4,1,0))
    if (bars){
        legend.DefaultArgs <- list(legend=names(object),col=col,cex=cex,bty="n",x="topleft")
        names.DefaultArgs <- list(cex=.7*par()$cex,y=c(-abs(diff(ylim))/15,-abs(diff(ylim))/25))
        frequencies.DefaultArgs <- list(cex=.7*par()$cex,percent=FALSE,offset=0)
    } else{
          legend.DefaultArgs <- list(legend=names(object),
                                     lwd=lwd,
                                     col=col,
                                     lty=lty,
                                     cex=cex,
                                     bty="n",
                                     y.intersp=1.3,
                                     x="topleft")
      }
    if(bars){
        if (type=="survival")
            legend.DefaultArgs$legend <- c("Predicted survival","Observed frequencies")
        else
            legend.DefaultArgs$legend <- c("Predicted risks","Observed frequencies")
    }
    lines.DefaultArgs <- list(type="l")
    abline.DefaultArgs <- list(lwd=1,col="red")
    if (missing(ylim)){
        if (showPseudo && !bars){
            ylim <- c(min(jack),max(jack))
        }
        else
            ylim <- c(0,1)
    }
    if (missing(xlim)){
        xlim <- c(0,1)
    }
    if (missing(xlab))
        if (bars)
            xlab <- ifelse(type=="survival","Survival groups","Risk groups")
        else
            xlab <- ifelse(type=="survival","Predicted survival probability","Predicted event probability")
    plot.DefaultArgs <- list(x=0,
                             y=0,
                             type = "n",
                             ylim = ylim,
                             xlim = xlim,
                             ylab=ylab,
                             xlab=xlab)
    barplot.DefaultArgs <- list(ylim = ylim,
                                col=col,
                                axes=FALSE,
                                ylab=ylab,
                                xlab=xlab,
                                beside=TRUE,
                                legend.text=NULL,
                                cex.axis=cex,
                                cex.lab=par()$cex.lab,
                                cex.names=cex)
    if (bars)
        control <- prodlim::SmartControl(call= list(...),
                                         keys=c("barplot","legend","axis2","abline","names","frequencies"),
                                         ignore=NULL,
                                         ignore.case=TRUE,
                                         defaults=list("barplot"=barplot.DefaultArgs,
                                             "abline"=abline.DefaultArgs,
                                             "legend"=legend.DefaultArgs,
                                             "names"=names.DefaultArgs,
                                             "frequencies"=frequencies.DefaultArgs,
                                             "axis2"=axis2.DefaultArgs),
                                         forced=list("abline"=list(h=0)),
                                         verbose=TRUE)
    else
        control <- prodlim::SmartControl(call= list(...),
                                         keys=c("plot","lines","legend","axis1","axis2"),
                                         ignore=NULL,
                                         ignore.case=TRUE,
                                         defaults=list("plot"=plot.DefaultArgs,
                                             "lines"=lines.DefaultArgs,
                                             "legend"=legend.DefaultArgs,
                                             "axis1"=axis1.DefaultArgs,
                                             "axis2"=axis2.DefaultArgs),
                                         forced=list("plot"=list(axes=FALSE),
                                             "axis1"=list(side=1)),
                                         verbose=TRUE)
    # }}}
    # {{{ splitmethod
    splitMethod <- resolvesplitMethod(splitMethod=splitMethod,B=B,N=NROW(data),M=M)
    k <- splitMethod$k
    B <- splitMethod$B
    N <- splitMethod$N
    NF <- length(object) 
    # }}}
    # {{{ ---------------------------Apparent predictions---------------------------
    apppred <- do.call("cbind",
                       lapply(1:NF,function(f){
                                  fit <- object[[f]]
                                  if (class(fit)[1] %in% c("numeric","double"))
                                      fit <- matrix(fit,ncol=1)
                                  apppred <- switch(model.type,
                                                    "competing.risks"={
                                                        p <- as.vector(do.call(predictHandlerFun,list(fit,newdata=data,times=time,cause=cause)))
                                                        if (class(fit)[[1]]%in% c("matrix","numeric")) p <- p[neworder]
                                                        p
                                                    },
                                                    "survival"={
                                                        p <- as.vector(do.call(predictHandlerFun,list(fit,newdata=data,times=time)))
                                                        if (class(fit)[[1]]%in% c("matrix","numeric")) 
                                                            p <- p[neworder]
                                                        p
                                                    },
                                                    "binary"={
                                                        p <- do.call(predictHandlerFun,list(fit,newdata=data))
                                                        if (class(fit)[[1]]%in% c("matrix","numeric")) p <- p[neworder]
                                                        p
                                                    })
                              }))
    colnames(apppred) <- names(object)
    if(pseudo==TRUE)
        apppred <- data.frame(jack=jack,apppred)
    else
        apppred <- data.frame(apppred)
    if (splitMethod$internal.name %in% c("noPlan")){
        predframe <- apppred
    }

    # }}}
    # {{{--------------k-fold and leave-one-out CrossValidation-----------------------
    if (splitMethod$internal.name %in% c("crossval","loocv")){
        groups <- splitMethod$index[,1,drop=TRUE]
        cv.list <- lapply(1:k,function(g){
                              if (verbose==TRUE) internalTalk(g,k)
                              id <- groups==g
                              train.k <- data[!id,,drop=FALSE]
                              val.k <- data[id,,drop=FALSE]
                              model.pred <- lapply(1:NF,function(f){
                                                       extraArgs <- giveToModel[[f]]
                                                       fit <- object[[f]]
                                                       fit.k <- internalReevalFit(object=fit,data=train.k,step=paste("CV group=",k),silent=FALSE,verbose=verbose)
                                                       switch(model.type,
                                                              "competing.risks"={do.call(predictHandlerFun,list(object=fit.k,newdata=val.k,times=time,cause=cause))},
                                                              "survival"={
                                                                  p <- do.call(predictHandlerFun,c(list(object=fit.k,newdata=val.k,times=time),extraArgs))
                                                                  p
                                                              },
                                                              "binary"={
                                                                  p <- do.call(predictHandlerFun,list(object=fit.k,newdata=val.k))
                                                                  p
                                                              })
                                                   })
                              model.pred
                          })
        predframe <- do.call("cbind",lapply(1:NF,function(f){
                                                pred <- do.call("rbind",lapply(cv.list,function(x)x[[f]]))
                                                if (splitMethod$internal.name!="loocv"){
                                                    pred <- pred[order(order(groups)),]
                                                }
                                                pred
                                            }))
        colnames(predframe) <- names(object)
        if(pseudo==TRUE)
            predframe <- cbind(data.frame(jack=jack),predframe)
        ## predframe <- na.omit(predframe)
    }
    # }}} 
    # {{{ ----------------------BootstrapCrossValidation----------------------
  
    if (splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632")){
        if (splitMethod$internal.name %in% c("Boot632plus","Boot632")){
            stop("Don't know how to do the 632(+) for the calibration curve.")
        }
        ResampleIndex <- splitMethod$index
        ## predframe <- do.call("rbind",lapply(1:B,function(b){
        ## predframe <- matrix
        pred.list <- parallel::mclapply(1:B,function(b){
                                            if (verbose==TRUE) internalTalk(b,B)
                                            jackRefit <- FALSE
                                            vindex.b <- match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0
                                            val.b <- data[vindex.b,,drop=FALSE]
                                            if (jackRefit){
                                                margFit.b <- prodlim::prodlim(margForm,data=val.b)
                                                jack.b <- prodlim::jackknife(margFit.b,cause=cause,times=time)
                                            }
                                            else{
                                                jack.b <- jack[match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0]
                                            }
                                            train.b <- data[ResampleIndex[,b],,drop=FALSE]
                                            frame.b <- data.frame(jack=jack.b)
                                            bootpred <- do.call("cbind",lapply(1:NF,function(f){
                                                                                   fit <- object[[f]]
                                                                                   fit.b <- internalReevalFit(object=fit,data=train.b,step=b,silent=FALSE,verbose=verbose)
                                                                                   extraArgs <- giveToModel[[f]]
                                                                                   try2predict <- try(pred.b <- switch(model.type,
                                                                                                                       "competing.risks"={do.call(predictHandlerFun,list(object=fit.b,newdata=val.b,times=time,cause=cause))},
                                                                                                                       "survival"={
                                                                                                                           p <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=val.b,times=time),extraArgs))
                                                                                                                           p
                                                                                                                       },
                                                                                                                       "binary"={
                                                                                                                           p <- do.call(predictHandlerFun,list(object=fit.b,newdata=val.b))
                                                                                                                           p
                                                                                                                       }),silent=TRUE)
                                                                                   if (inherits(try2predict,"try-error")==TRUE){
                                                                                       rep(NA,NROW(val.b))
                                                                                   }else{
                                                                                        pred.b
                                                                                    }
                                                                               }))
                                            colnames(bootpred) <- names(object)
                                            cbind(frame.b,bootpred)
                                        },mc.cores=cores)
        predframe <- do.call("rbind",pred.list)
        rm(pred.list)
    }
    # }}}
    # {{{ smoothing

    method <- match.arg(method,c("quantile","nne"))
    getXY <- function(f){
        if(pseudo==TRUE){
            p <- predframe[,f+1]
            jackF <- predframe[,1]
        }else{
             p <- predframe[,f]
         }
        switch(method,
               "quantile"={
                   if (length(q)==1)
                       groups <- quantile(p,seq(0,1,1/q))
                   else{
                       groups <- q
                   }
                   xgroups <- (groups[-(length(groups))]+groups[-1])/2
                   pcut <- cut(p,groups,include.lowest=TRUE)
                   if (pseudo==TRUE){
                       plotFrame=data.frame(Pred=tapply(p,pcut,mean),Obs=pmin(1,pmax(0,tapply(jackF,pcut,mean))))
                       attr(plotFrame,"quantiles") <- groups
                       plotFrame
                   }
                   else{
                       form.pcut <- update(formula,paste(".~pcut"))
                       if ("data.table" %in% class(data))
                           pdata <- cbind(data[,all.vars(update(formula,".~1")),drop=FALSE,with=FALSE],pcut=pcut)
                       else
                           pdata <- cbind(data[,all.vars(update(formula,".~1")),drop=FALSE],pcut=pcut)
                       y <- unlist(predict(f <- prodlim::prodlim(form.pcut,data=pdata),
                                           cause=cause,
                                           newdata=data.frame(pcut=levels(pcut)),
                                           times=time,
                                           type=ifelse(model.type=="survival","surv","cuminc")))
                       ## Is it ok to extrapolate into the future??
                       if (model.type=="survival")
                           y[is.na(y)] <- min(y,na.rm=TRUE)
                       else
                           y[is.na(y)] <- max(y,na.rm=TRUE)
                       plotFrame=data.frame(Pred=tapply(p,pcut,mean),Obs=y)
                       attr(plotFrame,"quantiles") <- groups
                       plotFrame
                   }
               },
               "nne"={
                   if (pseudo==TRUE){
                       ## Round probabilities to 2 digits
                       ## to avoid memory explosion ...
                       ## a difference in the 3 digit should
                       ## not play a role for the patient.
                       if (round==TRUE){
                           if (!is.null(bandwidth) && bandwidth>=1){
                               ## message("No need to round predicted probabilities to calculate calibration in the large")
                           } else{
                                 p <- round(p,2)
                             }
                       }
                       p <- na.omit(p)
                       if (no <- length(attr(p,"na.action")))
                           warning("calPlot: removed ",no," missing values in risk prediction.",call.=FALSE,immediate.=TRUE)
                       if (is.null(bandwidth)){
                           if (length(p)>length(apppred[,f+1])){
                               bw <- prodlim::neighborhood(apppred[,f+1])$bandwidth
                           }else{
                                bw <- prodlim::neighborhood(p)$bandwidth
                            }
                       } else{
                             bw <- bandwidth
                         }
                       if (bw>=1){
                           ## calibration in the large
                           plotFrame <- data.frame(Pred=mean(p),Obs=mean(jackF))
                       } else{
                             nbh <- prodlim::meanNeighbors(x=p,y=jackF,bandwidth=bw)
                             plotFrame <- data.frame(Pred=nbh$uniqueX,Obs=nbh$averageY)
                         }
                       attr(plotFrame,"bandwidth") <- bw
                       plotFrame
                   }else{
                        form.p <- update(formula,paste(".~p"))
                        if ("data.table" %in% class(data))
                            pdata <- cbind(data[,all.vars(update(formula,".~1")),drop=FALSE,with=FALSE],p=p)
                        else
                            pdata <- cbind(data[,all.vars(update(formula,".~1")),drop=FALSE],p=p)
                        y <- unlist(predict(prodlim::prodlim(form.p,data=pdata),
                                            cause=cause,
                                            newdata=data.frame(p=sort(p)),
                                            times=time,
                                            type=ifelse(type=="survival","surv","cuminc")))
                        plotFrame <- data.frame(Pred=sort(p),Obs=y)
                        plotFrame
                    }
               })
    }
    plotFrames <- lapply(1:NF,function(f){getXY(f)})
    names(plotFrames) <- names(object)
    # }}}
    # {{{ plot and/or invisibly output the results
    if (bars){
        if (model.type=="survival" && type=="risk")
            plotFrames[[1]] <- plotFrames[[1]][NROW(plotFrames[[1]]):1,]
        if ((is.logical(names[1]) && names[1]==TRUE)|| names[1] %in% c("quantiles.labels","quantiles")){
            qq <- attr(plotFrames[[1]],"quantiles")
            if (model.type=="survival" && type=="risk")
                qq <- rev(1-qq)
            if (names[1]=="quantiles.labels"){
                pp <- seq(0,1,1/q)
                names <- paste0("(",
                                sprintf("%1.0f",100*pp[-length(pp)]),",",
                                sprintf("%1.0f",100*pp[-1]),
                                ")\n",
                                sprintf("%1.1f",100*qq[-length(qq)])," - ",
                                sprintf("%1.1f",100*qq[-1]))
            }
            else 
                names <- paste0(sprintf("%1.1f",100*qq[-length(qq)])," - ",
                                sprintf("%1.1f",100*qq[-1]))
        }
    }
    summary <- list(n=NROW(data))
    if (model.type%in%c("survival","competing.risks"))
        summary <- c(summary,list("Event"=table(response[response[,"status"]!=0 & response[,"time"]<=time,ifelse(model.type=="survival","status","event")]),
                                  "Lost"=sum(response[,"status"]==0 & response[,"time"]<=time),
                                  "Event.free"=NROW(response[response[,"time"]>time,])))
    out <- list(plotFrames=plotFrames,
                predictions=apppred,
                time=time,
                cause=cause,
                pseudo=pseudo,
                summary=summary,
                control=control,
                legend=legend,
                bars=bars,
                diag=diag,
                add=add,
                legend=legend,
                names=names,
                method=method,
                model.type=model.type,
                type=type,
                axes=axes,
                percent=percent,
                hanging=hanging,
                showFrequencies=showFrequencies,
                col=col,
                ylim=ylim,
                xlim=xlim,
                ylab=ylab,
                xlab=xlab,
                lwd=lwd,
                lty=lty,
                pch=pch,
                lty=lty,
                NF=NF,
                pseudo.col=pseudo.col,
                pseudo.pch=pseudo.pch,
                showPseudo=showPseudo,
                jack.density=jack.density)
    if (method=="nne")
        out <- c(out,list(bandwidth=sapply(plotFrames,function(x)attr(x,"bandwidth"))))
    class(out) <- "calibrationPlot"
    if (plot){
        plot.calibrationPlot(out)
    }
    invisible(out)
    # }}}
}
