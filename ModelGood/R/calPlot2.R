##' Calibration plots for risk prediction models in for a binary endpoint
##'
##' For method "nne" the optimal bandwidth with respect to is obtained with
##' the function \code{\link{dpik}} from the package \code{KernSmooth} 
##' for a box kernel function. 
##' 
##' @title Calibration plots for binary data
##' @param object A named list of prediction models, where allowed
##' entries are (1) R-objects for which a \link{predictStatusProb}
##' method exists (see details), (2) a \code{call} that evaluates to
##' such an R-object (see examples), (3) a matrix with predicted
##' probabilities having as many rows as \code{data} in one column.
##' For cross-validation all objects in this list must include their
##' \code{call}.
##' @param formula A survival or event history formula. The left hand
##' side is used to compute the expected event status. If
##' \code{formula} is \code{missing}, try to extract a formula from
##' the first element in object.
##' @param data A data frame in which to validate the prediction
##' models and to fit the censoring model.  If \code{data} is missing,
##' try to extract a data set from the first element in object.
##' @param splitMethod Defines the internal validation design:
##'    
##'    \code{none/noPlan}: Assess the models in the give \code{data},
##'    usually either in the same data where they are
##'    fitted, or in independent test data. 
##'    
##'    \code{BootCv}: Bootstrap cross validation. The prediction models are
##'      trained on \code{B} bootstrap samples, that are either drawn with
##'      replacement of the same size as the original data or
##'      without replacement from \code{data} of the size \code{M}.
##'      The models are assessed in the observations that are NOT in the
##'      bootstrap sample.
##'      
##' @param B The number of cross-validation steps.
##' @param M The size of the subsamples for cross-validation.
##' @param showY If \code{TRUE} the observed data are shown as dots on
##' the plot.
##' @param method The method for estimating the calibration curve(s):
##'
##' \code{"nne"}: The expected event status is obtained in the nearest neighborhood around the predicted event probabilities.
##'
##' \code{"quantile"}: The expected event status is obtained in groups defined by quantiles of the predicted event probabilities.
##' 
##' @param round If \code{TRUE} predicted probabilities are rounded to
##' two digits before smoothing. This may have a considerable effect
##' on computing efficiency in large data sets.
##' @param bandwidth The bandwidth for \code{method="nne"}
##' @param q The number of quantiles for \code{method="quantile"}.
##' @param density Gray scale for observations.
##' @param add If \code{TRUE} the line(s) are added to an existing
##' plot.
##' @param diag If \code{FALSE} no diagonal line is drawn.
##' @param legend If \code{FALSE} no legend is drawn.
##' @param axes If \code{FALSE} no axes are drawn.
##' @param xlim Limits of x-axis.
##' @param ylim Limits of y-axis.
##' @param xlab Label for y-axis.
##' @param ylab Label for x-axis.
##' @param col Vector with colors, one for each element of
##' object. Passed to \code{\link{lines}}.
##' @param lwd Vector with line widths, one for each element of
##' object. Passed to \code{\link{lines}}.
##' @param lty lwd Vector with line style, one for each element of
##' object. Passed to \code{\link{lines}}.
##' @param pch Passed to \code{\link{points}}.
##' @param cause For competing risks models, the cause of failure or
##' event of interest
##' @param percent If TRUE axes labels are multiplied by 100 and thus
##' interpretable on a percent scale.
##' @param giveToModel List of with exactly one entry for each entry
##' in \code{object}. Each entry names parts of the value of the
##' fitted models that should be extracted and added to the value.
##' @param na.action Passed to \code{\link{model.frame}}
##' @param cores Number of cores for parallel computing.
##' Passed as the value of the argument \code{mc.cores}
##' when calling \code{\link{mclapply}}.
##' @param verbose if \code{TRUE} report details of the progress,
##' e.g. count the steps in cross-validation.
##' @param ... Used to control the subroutines: plot, axis, lines,
##' legend. See \code{\link{SmartControl}}.
##' @return list with elements: time, Frame and bandwidth (NULL for method quantile).
#' @references
#'  TA Gerds, PA Andersen, and Kattan MW. Calibration plots for risk prediction
#'  models in the presence of competing risks. Statistics in Medicine, page to
#'  appear, 2014.
#'  
##' @examples
#' set.seed(40)
#' N=40
#' Y=rbinom(N,1,.5)
#' X1=rnorm(N)
#' X1[Y==1]=rnorm(sum(Y==1),mean=rbinom(sum(Y==1),1,.5))
#' X2=rnorm(N)
#' X2[Y==0]=rnorm(sum(Y==0),mean=rbinom(sum(Y==0),3,.5))
#' dat <- data.frame(Y=Y,X1=X1,X2=X2)
#' lm1 <- glm(Y~X1,data=dat,family="binomial")
#' lm2 <- glm(Y~X2,data=dat,family="binomial")
#' calPlot2(list(lm1,lm2),data=dat)
##'
##'  
##' 
##' @author Thomas Alexander Gerds
##' @export
calPlot2 <- function(object,
                     formula,
                     data,
                     splitMethod="none",
                     B=1,
                     M,
                     showY,
                     method="nne",
                     round=TRUE,
                     bandwidth=NULL,
                     q=10,
                     density=55,
                     add=FALSE,
                     diag=!add,
                     legend=!add,
                     axes=!add,
                     xlim,
                     ylim,
                     xlab = "Predicted event probability",
                     ylab = "Observed proportion",
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
                     ...){

    if (missing(showY)) showY <- TRUE
    
  # {{{ find number of objects and lines
  cobj=class(object)[[1]]
  if (cobj!="list"){
    object <- list(object)
  }
  if (is.null(names(object))){
    names(object) <- sapply(object,function(o)class(o)[1])
  }
  else{
    names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
  }
  names(object) <- make.names(names(object),unique=TRUE)
  NF <- length(object)

    # }}}
    # {{{ lines types
    if (missing(lwd)) lwd <- rep(3,NF)
    if (missing(col)) col <- 1:NF
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
        formula <- eval(object[[1]]$call$formula)
        if (match("formula",class(formula),nomatch=0)==0)
            stop("Argument formula is missing.")
        else if (verbose)
            warning("Argument formula is missing. I use the formula from the call to the first model instead.")
    }
    m <- model.frame(formula,data,na.action=na.action)
    response <- model.response(m)
    stopifnot(length(unique(response))==2)
    predictHandlerFun <- "predictStatusProb"
    if (is.factor(response))
        jack <- as.numeric(response==levels(response)[2])
    else
        jack <- as.numeric(response)
    # }}}
    # {{{ call smartControls
    axis1.DefaultArgs <- list(side=1,las=1,at=seq(0,1,.25))
    axis2.DefaultArgs <- list(side=2,las=2,at=seq(0,1,.25),mgp=c(4,1,0))
    legend.DefaultArgs <- list(legend=names(object),lwd=lwd,col=col,lty=lty,cex=1.5,bty="n",y.intersp=1.3,x="topleft")
    lines.DefaultArgs <- list(type="l")
    if (missing(ylim)){
        if (showY){
            ylim <- c(min(jack),max(jack))
        }
        else
            ylim <- c(0,1)
    }
    if (missing(xlim)){
        xlim <- c(0,1)
    }
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,ylab="",xlab=xlab)
    smartA <- prodlim::SmartControl(call= list(...),keys=c("plot","lines","legend","axis1","axis2"),ignore=NULL,ignore.case=TRUE,defaults=list("plot"=plot.DefaultArgs,"lines"=lines.DefaultArgs,"legend"=legend.DefaultArgs,"axis1"=axis1.DefaultArgs,"axis2"=axis2.DefaultArgs),forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),verbose=TRUE)

    # }}}
    # {{{ splitmethod
    splitMethod <- MgSplitMethods(splitMethod=splitMethod,B=B,N=NROW(data),M=M)
    k <- splitMethod$k
    B <- splitMethod$B
    N <- splitMethod$N
    NF <- length(object)
    # }}}
    # {{{ cv, predictions and expectations
    # {{{ ---------------------------Apparent predictions---------------------------
    apppred <- do.call("cbind",lapply(1:NF,function(f){
        if (class(object[[f]][[1]])=="matrix"){
            apppred <- object[[f]][[1]]
        } else{
            apppred <- do.call(predictHandlerFun,list(object[[f]],newdata=data))
        }
    }))
    colnames(apppred) <- names(object)
    apppred <- data.frame(jack=jack,apppred)
    if (splitMethod$internal.name %in% c("noSplitMethod")){
        predframe <- apppred
    }
    # }}}
    # {{{--------------k-fold and leave-one-out CrossValidation-----------------------
    if (splitMethod$internal.name %in% c("crossval","loocv")){
        groups <- splitMethod$index[,1,drop=TRUE]
        cv.list <- lapply(1:k,function(g){
            if (verbose==TRUE) MgTalk(g,k)
            id <- groups==g
            train.k <- data[!id,,drop=FALSE]
            val.k <- data[id,,drop=FALSE]
            model.pred <- lapply(1:NF,function(f){
                extraArgs <- giveToModel[[f]]
                fit.k <- MgRefit(object=object[[f]],data=train.k,step=paste("CV group=",k),silent=FALSE,verbose=verbose)
                do.call(predictHandlerFun,list(object=fit.k,newdata=val.k))
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
            if (verbose==TRUE) MgTalk(b,B)
            jackRefit <- FALSE
            vindex.b <- match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0
            val.b <- data[vindex.b,,drop=FALSE]
            Y.b <- jack[match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0]
            train.b <- data[ResampleIndex[,b],,drop=FALSE]
            frame.b <- data.frame(jack=Y.b)
            bootpred <- do.call("cbind",lapply(1:NF,function(f){
                fit.b <- MgRefit(object=object[[f]],data=train.b,step=b,silent=FALSE,verbose=verbose)
                extraArgs <- giveToModel[[f]]
                try2predict <- try(pred.b <- do.call(predictHandlerFun,list(object=fit.b,newdata=val.b)),silent=TRUE)
                if (inherits(try2predict,"try-error")==TRUE){
                    rep(NA,NROW(val.b))
                }
                pred.b
            }))
            colnames(bootpred) <- names(object)
            cbind(frame.b,bootpred)
        },mc.cores=cores)
        predframe <- do.call("rbind",pred.list)
        rm(pred.list)
    }
  
    # }}}
    # }}}
    # {{{ smoothing
    method <- match.arg(method,c("quantile","nne"))
    outcome <- "pseudo"
    plotFrames <- lapply(1:NF,function(f){
        p <- predframe[,f+1]
        jackF <- predframe[,1]
        switch(method,
               "quantile"={
                   groups <- quantile(p,seq(0,1,1/q))
                   xgroups <- (groups[-(q+1)]+groups[-1])/2
                   plotFrame=data.frame(x=xgroups,y=tapply(jackF,cut(p,groups,include.lowest=TRUE),mean))
               },
               "nne"={
                   if (outcome=="pseudo"){
                       ## Round probabilities to 2 digits
                       ## to avoit memory explosion ...
                       ## a difference in the 3 digit should
                       ## not play a role for the patient.
                       if (round==TRUE) p <- round(p,2)
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
                       ## print(bw)
                       nbh <- prodlim::meanNeighbors(x=p,y=jackF,bandwidth=bw)
                       plotFrame <- data.frame(x=nbh$uniqueX,y=nbh$averageY)
                       attr(plotFrame,"bandwidth") <- bw
                       plotFrame
                   }
               })
    })
    # }}}
    # {{{ plot an empty frame
    if (add==FALSE){
        do.call("plot",smartA$plot)
        if (axes){
            if (percent){
                smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
                ## if (min(jack)<0) smartA$axis2$labels[1] <- ""
                ## if (max(jack)>1) smartA$axis2$labels[length(smartA$axis2$labels)] <- ""
                smartA$axis1$labels <- paste(100*smartA$axis1$at,"%")
            }
            do.call("axis",smartA$axis1)
            mgp2 <- smartA$axis2$mgp
            if (length(mgp2)>0){
                oldmgp <- par()$mgp
                par(mgp=mgp2)
                smartA$axis2 <- smartA$axis2[-match("mgp",names(smartA$axis2),nomatch=0)]
                title(ylab=ylab)
            }
            ## print(par()$mgp)
            do.call("axis",smartA$axis2)
            if (length(mgp2)>0){
                par(mgp=oldmgp)
            }
        }
    }
    if (diag){
        segments(x0=0,y0=0,x1=1,y1=1,col="gray77",lwd=2,xpd=FALSE)
    }
    ##   do.call("abline",c(list(a=0,b=1),list(col="gray77",lwd=2,xpd=FALSE)))
    # }}}
    # {{{ add lines and observations
    nix <- lapply(1:NF,function(f){
        plotFrame <- plotFrames[[f]]
        ## calibration in the large
        if (!is.null(bandwidth) && bandwidth>=1){
            with(na.omit(plotFrame),points(mean(x),mean(y),col=col[f],pch=16,cex=2))
        }else{
            with(na.omit(plotFrame),lines(x,y,col=col[f],lwd=lwd[f],lty=lty[f],type=ifelse(method=="quantile","b","l")))
        }
        ccrgb=as.list(col2rgb(col[f],alpha=TRUE))
        names(ccrgb) <- c("red","green","blue","alpha")
        ccrgb$alpha <- density
        Y.col <- do.call("rgb",c(ccrgb,list(max=255)))
        if (showY) {
            points(apppred[,f+1],apppred[,1],col=Y.col)
        }
    })

  # }}}
  # {{{ legend
  ## if (missing(legend)) legend=ifelse(length(object)==1,FALSE,TRUE)
  ## if (missing(legend.legend)) legend.legend=names(object)
  if(legend){
    do.call("legend",smartA$legend)
  }
  ## if (legend)
  ## legend(0,1,legend=legend.legend,lwd=lwd,col=col,bty="n")

  # }}}
  # {{{ invisibly output the jackknife pseudo-values
    out <- list(Frame=predframe,bandwidth=sapply(plotFrames,function(x)attr(x,"bandwidth")))
  invisible(out)
  # }}}
}
