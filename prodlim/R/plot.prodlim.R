# {{{ Header
#' Plotting event probabilities over time
#' 
#' Function to plot survival and cumulative incidence curves against time.
#' 
#' From version 1.1.3 on the arguments legend.args, atrisk.args, confint.args
#' are obsolete and only available for backward compatibility. Instead
#' arguments for the invoked functions \code{atRisk}, \code{legend},
#' \code{confInt}, \code{markTime}, \code{axis} are simply specified as
#' \code{atrisk.cex=2}. The specification is not case sensitive, thus
#' \code{atRisk.cex=2} or \code{atRISK.cex=2} will have the same effect.  The
#' function \code{axis} is called twice, and arguments of the form
#' \code{axis1.labels}, \code{axis1.at} are used for the time axis whereas
#' \code{axis2.pos}, \code{axis1.labels}, etc. are used for the y-axis.
#' 
#' These arguments are processed via \code{\dots{}} of \code{plot.prodlim} and
#' inside by using the function \code{SmartControl}.  Documentation of these
#' arguments can be found in the help pages of the corresponding functions.
#' 
#' @aliases plot.prodlim lines.prodlim
#' @param x an object of class `prodlim' as returned by the
#' \code{prodlim} function.
#' @param type Either \code{"surv"} or \code{"cuminc"} controls what part of the object is plotted. 
# Defaults to \code{object$type}.  
#' @param cause determines the cause of the cumulative incidence
#' function.  Currently one cause is allowed at a time, but you may
#' call the function again with add=TRUE to add the lines of the other
#' causes.
#' @param select Select which lines to plot. This can be used when there
#' are many strata or many competing risks to select a subset of the lines.
#' However, a more clean way to select covariate stratat is to use argument \code{newdata}.
#' Another application is when there are many competing risks and it is desired (for the stacked plot)
#'  to stack and show only a subset of the cumulative incidence functions. 
#' @param newdata a data frame containing covariate strata for which to show 
#' curves. When omitted element \code{X} of object
#' \code{x} is used.
#' @param add if \code{TRUE} curves are added to an existing plot.
#' @param col color for curves. Default is \code{1:number(curves)}
#' @param lty line type for curves. Default is 1.
#' @param lwd line width for all curves. Default is 3.
#' @param ylim limits of the y-axis
#' @param xlim limits of the x-axis
#' @param xlab label for the x-axis
#' @param ylab label for the y-axis
#' @param legend if TRUE a legend is plotted by calling the function
#' legend.  Optional arguments of the function \code{legend} can be
#' given in the form \code{legend.x=val} where x is the name of the
#' argument and val the desired value. See also Details.
#' @param logrank If TRUE, the logrank p-value will be extracted from
#' a call to \code{survdiff} and added to the legend. This works only
#' for survival models, i.e. Kaplan-Meier with discrete predictors.
#' @param marktime if TRUE the curves are tick-marked at right
#' censoring times by invoking the function \code{markTime}. Optional
#' arguments of the function \code{markTime} can be given in the form
#' \code{confint.x=val} as with legend. See also Details.
#' @param confint if TRUE pointwise confidence intervals are plotted
#' by invoking the function \code{confInt}. Optional arguments of the
#' function \code{confInt} can be given in the form
#' \code{confint.x=val} as with legend.  See also Details.
#' @param automar If TRUE the function trys to find suitable values
#' for the figure margins around the main plotting region.
#' @param atrisk if TRUE display numbers of subjects at risk by
#' invoking the function \code{atRisk}. Optional arguments of the
#' function \code{atRisk} can be given in the form \code{atrisk.x=val}
#' as with legend. See also Details.
#' @param timeOrigin Start of the time axis
#' @param axes If true axes are drawn. See details.
#' @param background If \code{TRUE} the background color and grid
#' color can be controlled using smart arguments SmartControl, such as
#' background.bg="yellow" or background.bg=c("gray66","gray88").  The
#' following defaults are passed to \code{background} by
#' \code{plot.prodlim}: horizontal=seq(0,1,.25), vertical=NULL,
#' bg="gray77", fg="white".  See \code{background} for all arguments,
#' and the examples below.
#' @param percent If true the y-axis is labeled in percent.
#' @param minAtrisk Integer. Show the curve only until the number
#' at-risk is at least \code{minAtrisk}
#' @param limit When newdata is not specified and the number of lines
#' in element \code{X} of object \code{x} exceeds limits, only the
#' results for covariate constellations of the first, the middle and
#' the last row in \code{X} are shown. Otherwise all lines of \code{X}
#' are shown.
#' @param ... Parameters that are filtered by
#' \code{\link{SmartControl}} and then passed to the functions
#' \code{\link{plot}}, \code{\link{legend}}, \code{\link{axis}},
#' \code{\link{atRisk}}, \code{\link{confInt}},
#' \code{\link{markTime}}, \code{\link{backGround}}
#' @return The (invisible) object.
#' @note Similar functionality is provided by the function
#' \code{\link{plot.survfit}} of the survival library
#' @author Thomas Alexander Gerds <tag@@biostat.ku.dk>
#' @seealso \code{\link{plot}}, \code{\link{legend}}, \code{\link{axis}},
#' \code{\link{prodlim}},\code{\link{plot.Hist}},\code{\link{summary.prodlim}},
#' \code{\link{neighborhood}}, \code{\link{atRisk}}, \code{\link{confInt}},
#' \code{\link{markTime}}, \code{\link{backGround}}
#' @keywords survival
##' @examples
##' ## simulate right censored data from a two state model 
##' set.seed(100)
##' dat <- SimSurv(100)
##' # with(dat,plot(Hist(time,status)))
##' 
##' ### marginal Kaplan-Meier estimator
##' kmfit <- prodlim(Hist(time, status) ~ 1, data = dat)
##' plot(kmfit)
##' 
##' # change time range
##' plot(kmfit,xlim=c(0,4))
##' 
##' # change scale of y-axis
##' plot(kmfit,percent=FALSE)
##' 
##' # mortality instead of survival
##' plot(kmfit,type="cuminc")
##' 
##' # change axis label and position of ticks
##' plot(kmfit,
##'      xlim=c(0,10),
##'      axis1.at=seq(0,10,1),
##'      axis1.labels=0:10,
##'      xlab="Years",
##'      axis2.las=2,
##'      atrisk.at=seq(0,10,2.5),
##'      atrisk.title="")
##' 
##' # change background color
##' plot(kmfit,
##'      xlim=c(0,10),
##'      confint.citype="shadow",
##'      col=1,
##'      axis1.at=0:10,
##'      axis1.labels=0:10,
##'      xlab="Years",
##'      axis2.las=2,
##'      atrisk.at=seq(0,10,2.5),
##'      atrisk.title="",
##'      background=TRUE,
##'      background.fg="white",
##'      background.horizontal=seq(0,1,.25/2),
##'      background.vertical=seq(0,10,2.5),
##'      background.bg=c("gray88"))
##' 
##' # change type of confidence limits
##' plot(kmfit,
##'      xlim=c(0,10),
##'      confint.citype="dots",
##'      col=4,
##'      background=TRUE,
##'      background.bg=c("white","gray88"),
##'      background.fg="gray77",
##'      background.horizontal=seq(0,1,.25/2),
##'      background.vertical=seq(0,10,2))
##' 
##' 
##' ### Kaplan-Meier in discrete strata
##' kmfitX <- prodlim(Hist(time, status) ~ X1, data = dat)
##' plot(kmfitX)
##' # move legend
##' plot(kmfitX,legend.x="bottomleft",atRisk.cex=1.3,
##'      atrisk.title="No. subjects")
##' 
##' ## Control the order of strata
##' ## since version 1.5.1 prodlim does  obey the order of
##' ## factor levels
##' dat$group <- factor(cut(dat$X2,c(-Inf,0,0.5,Inf)),
##'                     labels=c("High","Intermediate","Low"))
##' kmfitG <- prodlim(Hist(time, status) ~ group, data = dat)
##' plot(kmfitG)
##' 
##' ## relevel 
##' dat$group2 <- factor(cut(dat$X2,c(-Inf,0,0.5,Inf)),
##'                      levels=c("(0.5, Inf]","(0,0.5]","(-Inf,0]"),
##'                      labels=c("Low","Intermediate","High"))
##' kmfitG2 <- prodlim(Hist(time, status) ~ group2, data = dat)
##' plot(kmfitG2)
##' 
##' # add log-rank test to legend
##' plot(kmfitX,
##'      atRisk.cex=1.3,
##'      logrank=TRUE,
##'      legend.x="topright",
##'      atrisk.title="at-risk")
##' 
##' # change atrisk labels
##' plot(kmfitX,
##'      legend.x="bottomleft",
##'      atrisk.title="Patients",
##'      atrisk.cex=0.9,
##'      atrisk.labels=c("X1=0","X1=1"))
##' 
##' # multiple categorical factors
##' 
##' kmfitXG <- prodlim(Hist(time,status)~X1+group2,data=dat)
##' plot(kmfitXG,select=1:2)
##' 
##' ### Kaplan-Meier in continuous strata
##' kmfitX2 <- prodlim(Hist(time, status) ~ X2, data = dat)
##' plot(kmfitX2,xlim=c(0,10))
##' 
##' # specify values of X2 for which to show the curves 
##' plot(kmfitX2,xlim=c(0,10),newdata=data.frame(X2=c(-1.8,0,1.2)))
##' 
##' ### Cluster-correlated data
##' library(survival)
##' cdat <- cbind(SimSurv(20),patnr=sample(1:5,size=20,replace=TRUE))
##' kmfitC <- prodlim(Hist(time, status) ~ cluster(patnr), data = cdat)
##' plot(kmfitC)
##' plot(kmfitC,atrisk.labels=c("Units","Patients"))
##' 
##' kmfitC2 <- prodlim(Hist(time, status) ~ X1+cluster(patnr), data = cdat)
##' plot(kmfitC2)
##' plot(kmfitC2,atrisk.labels=c("Teeth","Patients","Teeth","Patients"),
##'      atrisk.col=c(1,1,2,2))
##' 
##' 
##' ### Cluster-correlated data with strata
##' n = 50
##' foo = runif(n)
##' bar = rexp(n)
##' baz = rexp(n,1/2)
##' d = stack(data.frame(foo,bar,baz))
##' d$cl = sample(10, 3*n, replace=TRUE)
##' fit = prodlim(Surv(values) ~ ind + cluster(cl), data=d)
##' plot(fit)
##' 
##' 
##' ## simulate right censored data from a competing risk model 
##' datCR <- SimCompRisk(100)
##' with(datCR,plot(Hist(time,event)))
##' 
##' ### marginal Aalen-Johansen estimator
##' ajfit <- prodlim(Hist(time, event) ~ 1, data = datCR)
##' plot(ajfit) # same as plot(ajfit,cause=1)
##' 
##' # cause 2
##' plot(ajfit,cause=2)
##' 
##' # both in one
##' plot(ajfit,cause=1)
##' plot(ajfit,cause=2,add=TRUE,col=2)
##' 
##' ### stacked plot
##' 
##' plot(ajfit,cause="stacked",select=2)
##' 
##' ### stratified Aalen-Johansen estimator
##' ajfitX1 <- prodlim(Hist(time, event) ~ X1, data = datCR)
##' plot(ajfitX1)
##' 
##' ## add total number at-risk to a stratified curve
##' ttt = 1:10
##' plot(ajfitX1,atrisk.at=ttt,col=2:3)
##' plot(ajfit,add=TRUE,col=1)
##' atRisk(ajfit,newdata=datCR,col=1,times=ttt,line=3,labels="Total")
##' 
##'
##' ## stratified Aalen-Johansen estimator in nearest neighborhoods
##' ## of a continuous variable
##' ajfitX <- prodlim(Hist(time, event) ~ X1+X2, data = datCR)
##' plot(ajfitX,newdata=data.frame(X1=c(1,1,0),X2=c(4,10,10)))
##' plot(ajfitX,newdata=data.frame(X1=c(1,1,0),X2=c(4,10,10)),cause=2)
##' 
##' ## stacked plot
##' 
##' plot(ajfitX,
##'      newdata=data.frame(X1=0,X2=0.1),
##'      cause="stacked",
##'      legend.title="X1=0,X2=0.1",
##'      legend.legend=paste("cause:",getStates(ajfitX$model.response)),
##'      plot.main="Subject specific stacked plot")
##'  
#' @export 
plot.prodlim <- function(x,
                         type,
                         cause=1,
                         select,
                         newdata,
                         add = FALSE,
                         col,
                         lty,
                         lwd,
                         ylim,
                         xlim,
                         xlab="Time",
                         ylab,
                         legend=TRUE,
                         logrank=FALSE,
                         marktime=FALSE,
                         confint=TRUE,
                         automar,
                         atrisk=ifelse(add,FALSE,TRUE),
                         timeOrigin=0,
                         axes=TRUE,
                         background=TRUE,
                         percent=TRUE,
                         minAtrisk=0,
                         limit=10,
                         ## time.scale="same",
                         ...){

    # }}}
    # {{{  backward compatibility
    ##   args=match.call(expand=TRUE)
    ##   args[[1]]=list
    allArgs <- match.call()
    if (missing(type)){
        type=allArgs[[match("what",names(allArgs))]]
    }
    # }}}
    # {{{  extracting a list of lines to draw

  cens.type <- x$cens.type    # uncensored, right or interval censored
  if (cens.type=="intervalCensored") {
      confint <- FALSE
      atrisk <- FALSE
  }
  model <- x$model                 # survival, competing risks or multi-state
  clusterp <- !is.null(x$clustervar)
  if (missing(type)||is.null(type)){
      type <- x$type
      ## type <- switch(model,"survival"="surv","competing.risks"="cuminc","multi.states"="hazard")
      ## if (!is.null(x$reverse) && x$reverse==TRUE && model=="survival") type <- "cuminc"
  }
  else
      type <- match.arg(type,c("surv","cuminc","hazard"))
  if (model=="competing.risks" && type=="surv")
      stop("To plot the event-free survival curve, please fit a suitable model: prodlim(Hist(time,status!=0)~....")
  
  if (cens.type=="intervalCensored")
      plot.times <- sort(unique(x$time[2,]))
  else{
      plot.times <- sort(unique(x$time))
      if (plot.times[1]>timeOrigin) plot.times <- c(timeOrigin,plot.times)
      else plot.times <- c(timeOrigin,plot.times[plot.times>timeOrigin])
  }
  if (length(x$clustervar)>0)
      nRisk <- x$n.risk[,1]
  else
      nRisk <- x$n.risk
  if (minAtrisk>0 && any(nRisk<=minAtrisk)){
      if (all(nRisk<=minAtrisk)){
          return(plot(0,0,type="n",xlim=c(min(plot.times), max(plot.times)),ylim=c(0, 1),axes=FALSE))
      }
      criticalTime <- min(x$time[nRisk<=minAtrisk])
      plot.times <- plot.times[plot.times<criticalTime]
      ## if (sum(nEvent[nRisk>minAtrisk])<=1)
  }
  if (missing(newdata)) {
      newdata <- x$X
      if (NROW(newdata)>limit)
          newdata <- newdata[c(1,round(median(1:NROW(newdata))),NROW(newdata)),,drop=FALSE]          
  }
  ## restrict plot.times to xlim
  if (!missing(xlim)){
      if (xlim[1]>plot.times[1]) plot.times <- plot.times[plot.times>=xlim[1]]
      if (xlim[2]<plot.times[length(plot.times)]) plot.times <- plot.times[plot.times<=xlim[2]]
  }
  ## if (missing(newdata) && NROW(newdata)>limit)
  ## newdata <- newdata[c(1,round(median(1:NROW(newdata))),NROW(newdata)),,drop=FALSE]
  ## browser()
  stacked <- cause[1]=="stacked"
  if (stacked){
      confint <- FALSE
      if (model!="competing.risks") stop("Stacked plot works only for competing risks models.")
      if (NROW(newdata)>1) stop("Stacked plot works only for one covariate stratum.")
  }else{
       if (length(cause)!=1){
           warning("Currently only the cumulative incidence of a single cause can be plotted in one go. Use argument add=TRUE to add the lines of the other causes. For now I use the first cause")
           cause <- cause[1]
       }
   }
  ## Y <- predict(x,times=plot.times,newdata=newdata,level.chaos=1,type=type,cause=cause,mode="list")
  startValue=ifelse(type=="surv",1,0)
  if (type=="hazard" && model!="survival")
      stats=list(c("cause.hazard",0))
  else
      stats=list(c(type,startValue))
  if (model=="survival" && type=="cuminc") {
      startValue=1
      stats=list(c("surv",startValue))
  }
  if (confint==TRUE)
      stats=c(stats,list(c("lower",startValue),c("upper",startValue)))
  if (x$cens.type=="intervalCensored"){
      stop("FIXME: There is no plot method implemented for intervalCensored data.")
  }
  if  (model=="competing.risks"){
      if (stacked) ## all causes
          cause <- attributes(x$model.response)$states
      else
          cause <- checkCauses(cause,x)
      sumX <- lifeTab(x,
                      times=plot.times,
                      cause=cause,
                      newdata=newdata,
                      stats=stats,
                      percent=FALSE)
  }
  else{
      sumX <- lifeTab(x,
                      times=plot.times,
                      newdata=newdata,
                      stats=stats,
                      percent=FALSE)
  }
  if (model=="competing.risks"){
      if (stacked == FALSE){
          sumX <- sumX[[cause]]
      } else {
            ## there is at most one stratum for each cause
            if (!is.null(newdata))
                sumX <- lapply(sumX,function(cc)cc[[1]])
        }
  }
  ## cover both no covariate and single newdata:
  if (!is.null(dim(sumX))) sumX <- list(sumX)
  if (model=="survival" && type=="cuminc"){
      Y <- lapply(sumX,function(x)1-x[,"surv"])
      names(Y) <- names(sumX)
      nlines <- length(Y)
  } else{
        Y <- lapply(sumX,function(x)x[,type])
        names(Y) <- names(sumX)
        if (!missing(select)){
            if (length(select)==1)
                Y <- Y[select]
            else
                Y <- Y[select]
        }
        nlines <- length(Y)
    }
  
  # }}}
  # {{{  getting default arguments for plot, atrisk, axes, legend, confint, marktime
  
  if (missing(ylab)) ylab <- switch(type,"surv"=ifelse(x$reverse==TRUE,"Censoring probability","Survival probability"),"cuminc"="Cumulative incidence","hazard"="Cumulative hazard")
  if (missing(xlab)) xlab <- "Time"
  if (missing(xlim)) xlim <- c(min(plot.times), max(plot.times))
  if (missing(ylim)) ylim <- c(0, 1)
  if (missing(lwd)) lwd <- rep(3,nlines)
  if (missing(col)) col <- 1:nlines
  if (missing(lty)) lty <- rep(1, nlines)
  if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  
  background.DefaultArgs <- list(xlim=xlim,
                                 ylim=ylim,
                                 horizontal=seq(ylim[1],ylim[2],diff(ylim)/4),
                                 vertical=NULL,
                                 bg="white",
                                 fg="gray88")
  axis1.DefaultArgs <- list()
  axis2.DefaultArgs <- list(at=seq(ylim[1],ylim[2],ylim[2]/4),side=2)
  lines.DefaultArgs <- list(type="s")
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = xlab,ylab = ylab)
  marktime.DefaultArgs <- list(x=Y,nlost=lapply(sumX,function(x)x[,"n.lost"]),times=plot.times,pch="I",col=col)
  if (length(Y)==1 && length(x$clustervar)==0){
      atriskDefaultLabels <- "Subjects: "
      atriskDefaultTitle <- ""
  }
  else{
      if (length(x$clustervar)>0){
          atriskDefaultTitle <- ""
          atriskDefaultLabels <- rep(paste(c("Subjects","Clusters"),": ",sep=""),
                                     nlines)
      }
      else{
          ## print(names(Y))
          if (model=="competing.risks" && stacked==TRUE){
              atriskDefaultTitle <- ""
              atriskDefaultLabels <- "Subjects: "
          }
          else{

              if ((length(grep("=",names(Y)))==length(names(Y)))){
                  atriskDefaultLabels <- paste(gsub("[ \t]*$","",sapply(strsplit(names(Y),"="),function(x)x[[2]])),
                                               ": ", sep="")
                  atriskDefaultTitle <- unique(sapply(strsplit(names(Y),"="),function(x)x[[1]]))
              }else{
                   atriskDefaultTitle <- ""
                   atriskDefaultLabels <- paste(gsub("[ \t]*$","",names(Y)),": ",sep="")
               }
          }
      }
      ## atriskDefaultLabels <- format(atriskDefaultLabels,justify="left")
      ## atriskDefaultTitle <- ""
  }
  atrisk.DefaultArgs <- list(x=x,
                             newdata=newdata,
                             interspace=1,
                             dist=.3,
                             col=col,
                             labelcol=1,
                             titlecol=1,
                             title=atriskDefaultTitle,
                             labels=atriskDefaultLabels,
                             times=seq(min(plot.times),max(plot.times),(max(plot.times)-min(plot.times))/10))
                             ## times=seq(0,min(x$maxtime,xlim[2]),min(x$maxtime,xlim[2])/10))
  if (!missing(select) && (!(model=="competing.risks" && stacked))){
      atrisk.DefaultArgs$newdata <- atrisk.DefaultArgs$newdata[select,,drop=FALSE]
  }
  legend.DefaultArgs <- list(legend=names(Y),
                             lwd=lwd,
                             col=col,
                             lty=lty,
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             trimnames=!match("legend.legend",names(allArgs),nomatch=0),
                             x="topright")
  if (stacked) {
      legend.DefaultArgs$title <- "Competing risks"
      legend.DefaultArgs$x <- "topleft"
  }

  if (NCOL(newdata)>1) legend.DefaultArgs$trimnames <- FALSE
  confint.DefaultArgs <- list(x=x,
                              newdata=newdata,
                              type=type,
                              citype="shadow",
                              times=plot.times,
                              cause=cause,
                              density=55,
                              col=col[1:nlines],
                              lwd=rep(2,nlines),
                              lty=rep(3,nlines))

  # }}}
# {{{  backward compatibility

  if (match("legend.args",names(allArgs),nomatch=FALSE)){
      legend.DefaultArgs <- c(args[[match("legend.args",names(allArgs),nomatch=FALSE)]],legend.DefaultArgs)
      legend.DefaultArgs <- legend.DefaultArgs[!duplicated(names(legend.DefaultArgs))]
  }
  if (match("confint.args",names(allArgs),nomatch=FALSE)){
      confint.DefaultArgs <- c(args[[match("confint.args",names(allArgs),nomatch=FALSE)]],confint.DefaultArgs)
      confint.DefaultArgs <- confint.DefaultArgs[!duplicated(names(confint.DefaultArgs))]
  }
  if (match("atrisk.args",names(allArgs),nomatch=FALSE)){
      atrisk.DefaultArgs <- c(args[[match("atrisk.args",names(allArgs),nomatch=FALSE)]],atrisk.DefaultArgs)
      atrisk.DefaultArgs <- atrisk.DefaultArgs[!duplicated(names(atrisk.DefaultArgs))]
  }
  if (length(list(...)) && match("legend.legend",names(list(...)),nomatch=FALSE) && any(sapply(newdata,is.factor))){
      message("Since version 1.5.1 prodlim obeys the order of factor levels.\nThis may break old code which explicitly defines the legend labels.")
  }
  
  smartA <- SmartControl(call=  list(...),
                         keys=c("plot","lines","atrisk","legend","confint","background","marktime","axis1","axis2"),
                         ignore=c("x","type","cause","newdata","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","marktime","confint","automar","atrisk","timeOrigin","percent","axes","atrisk.args","confint.args","legend.args"),
                         defaults=list("plot"=plot.DefaultArgs,"atrisk"=atrisk.DefaultArgs,"lines"=lines.DefaultArgs,"legend"=legend.DefaultArgs,"confint"=confint.DefaultArgs,"marktime"=marktime.DefaultArgs,"background"=background.DefaultArgs,"axis1"=axis1.DefaultArgs,"axis2"=axis2.DefaultArgs),
                         forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),
                         ignore.case=TRUE,
                         replaceDefaults=FALSE,
                         verbose=TRUE)

  # }}}
  # {{{  setting margin parameters
  if (atrisk==TRUE){
      oldmar <- par()$mar
      if (missing(automar) || automar==TRUE){
          ##        bottomMargin =  margin line (in 'mex' units) for xlab
          ##                        + distance of xlab from xaxis
          ##                        + distance of atrisk numbers from xlab
          ##                        + number of atrisk lines
          ##                        + one extra line below the bottom number atrisk line
          ##      leftSideMargin =  margin line + atrisk.lab
          bottomMargin <- par()$mgp[2] + smartA$atrisk$dist+ ifelse(clusterp,2,1)*nlines + 1
          ## smartA$atrisk$labels
          maxlabellen <- max(strwidth(c(smartA$atrisk$labels,smartA$atrisk$title),
                                      cex=smartA$atrisk$cex,
                                      units="inches"))
          maxlabellen <- pmax(maxlabellen * (par("mar")[2] / par("mai")[2]),par("mar")[2])
          leftMargin <- maxlabellen+2-par("mar")[2]
          newmar <- par()$mar + c(bottomMargin,leftMargin,0,0)
          par(mar=newmar)
      }
  }

  # }}}
# {{{  plot and backGround
  if (!add) {
    do.call("plot",smartA$plot)
    ##     if (background==TRUE && match("bg",names(smartA$background),nomatch=FALSE)){
    ## par(bg=smartA$background$bg)
    ##     }
    if (background==TRUE){
      do.call("backGround",smartA$background)
    }
  }
  # }}}
# {{{  axes

  if (!add) {
    if (axes){
      do.call("axis",smartA$axis1)
      if (percent & is.null(smartA$axis2$labels))
        smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
      do.call("axis",smartA$axis2)
    }
  }
  if (atrisk==TRUE) par(mar=oldmar) ## reset

  # }}}
  # {{{  pointwise confidence intervals
  if (confint==TRUE) {
      ## if (verbose==TRUE){print(smartA$confint)}
      do.call("confInt",smartA$confint)
  }
  # }}}
  # {{{  adding the lines 
  lines.type <- smartA$lines$type
  if (stacked==TRUE){
      if (length(Y)>1){
          nY <- names(Y)
          Y <- apply(do.call("rbind",Y),2,cumsum)
          Y <- lapply(1:nlines,function(i)Y[i,])
          names(Y) <- nY
      }
      ## names(Y) <- attr(x$model.response,"states")
      nix <- lapply(1:nlines, function(s) {
                        yyy <- Y[[s]]
                        ppp <- plot.times
                        pos.na <- is.na(yyy)
                        ppp <- ppp[!pos.na]
                        yyy <- yyy[!pos.na]
                        lines(x = ppp,y = yyy,type = lines.type,col = col[s],lty = lty[s],lwd = lwd[s])
                        cc <- dimColor(col[s],density=55)
                        ttt <- ppp
                        nt <- length(ttt)
                        ttt <- c(ttt,ttt)
                        uuu <- c(0,yyy[-nt],yyy)
                        if (s==1)
                            lll <- rep(0,nt*2)
                        else
                            lll <- c(0,Y[[s-1]][!pos.na][-nt],Y[[s-1]][!pos.na])
                        neworder <- order(ttt)
                        uuu <- uuu[neworder]
                        lll <- lll[neworder]
                        ttt <- sort(ttt)
                        polygon(x=c(ttt,rev(ttt)),y=c(lll,rev(uuu)),col=cc,border=NA)
                    })
  }else{
       nix <- lapply(1:nlines, function(s) {
                         lines(x = plot.times,
                               y = Y[[s]],
                               type = lines.type,
                               col = col[s],
                               lty = lty[s],
                               lwd = lwd[s])
                     })
   }
  # }}}
  # {{{  marks at the censored times

  if (marktime==TRUE){
    if (model %in% c("survival","competing.risks")){
      do.call("markTime",smartA$marktime)
    }
    else{
      message("Marking the curves at censored times is not yet available for multi-state models.")
    }
  }

# }}}
# {{{  adding the no. of individuals at risk

  if (atrisk==TRUE && !add){
      if (hit <- match("at",names(smartA$atrisk),nomatch=FALSE)){
          if (match("atrisk.times",names(list(...)),nomatch=FALSE)){
              warning("Atrisk argument clash: remove either 'atrisk.at' or 'atrisk.times'.")
          }
          else{
              names(smartA$atrisk)[hit] <- "times"
              smartA$atrisk <- smartA$atrisk[!duplicated(names(smartA$atrisk))]
          }
      }
      do.call("atRisk",smartA$atrisk)
  }
  # }}}
  # {{{  legend
  if(legend==TRUE && !add && !is.null(names(Y))){
      if (smartA$legend$trimnames==TRUE && (length(grep("=",smartA$legend$legend))==length(smartA$legend$legend))){
          smartA$legend$legend <- sapply(strsplit(smartA$legend$legend,"="),function(x)x[[2]])
          if (is.null(smartA$legend$title))
              smartA$legend$title <- unique(sapply(strsplit(names(Y),"="),function(x)x[[1]]))
      }
      smartA$legend <- smartA$legend[-match("trimnames",names(smartA$legend))]
      save.xpd <- par()$xpd
      if (logrank && model=="survival" && length(smartA$legend$legend)>1){
          ## formula.names <- try(all.names(formula),silent=TRUE)
          lrform <- x$call$formula
          if (lrform[[2]][[1]]==as.name("Hist"))
              lrform[[2]][[1]] <- as.name("Surv")
          ## require(survival)
          lrtest <- survival::survdiff(eval(lrform),data=eval(x$call$data))
          if (length(lrtest$n) == 1) {
              p <- 1 - pchisq(lrtest$chisq, 1)
          } else{
              if (is.matrix(x$obs)) {
                  etmp <- apply(lrtest$exp, 1, sum)
              }
              else {
                  etmp <- lrtest$exp
              }
              df <- (sum(1 * (etmp > 0))) - 1
              p <- 1 - pchisq(lrtest$chisq, df)
          }
          if (length(smartA$legend$title))
              smartA$legend$title <- paste(smartA$legend$title," Log-rank: p=",format.pval(p,digits=logrank,eps=0.0001))
          else
              smartA$legend$title <- paste(" Log-rank: ",format.pval(p,digits=logrank,eps=0.0001))
      }
      par(xpd=TRUE)
      do.call("legend",smartA$legend)
      par(xpd=save.xpd)
  }

# }}}
invisible(x)
}
