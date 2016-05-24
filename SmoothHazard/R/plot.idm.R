#' Plot method for an illness-death model 
#' 
#' Plot estimated baseline transition intensities from an object of class
#' \code{idm} optionally with confidence limits.
#' 
#' 
#' @param x a \code{idmWeib} class object (output from calling
#' \code{idm} with the (default) option \code{intensities}="Weib".
#' @param conf.int If TRUE show confidence limits
#' @param citype Type of confidence limits, can be "shadow" or "bars"
#' @param add If TRUE add to existing plot
#' @param axes If TRUE axes are drawn
#' @param col Color of the lines
#' @param lwd Width of the lines
#' @param lty Type of the lines
#' @param xlim Limits for x-axis
#' @param ylim Limits for y-axis
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param legend If TRUE a legend is drawn, which can be further controlled via \code{\link{SmartControl}}.
#' @param transition Choose one of the transition intensities: \code{c("01","02","12")}.
#' @param ... Passed to \code{\link{SmartControl}}
#' @return Print a plot of the baseline transition intensities of an
#' illness-death model estimated using a Weibull approach.
#' @seealso
#' \code{\link{print.idm}},\code{\link{summary.idm}},\code{\link{idm}},
#' @seealso \code{\link{idm}}
#' @keywords methods
#' @examples
#'
#' library(lava)
#' library(prodlim)
#' m <- idmModel(scale.lifetime=1/10,scale.illtime=1/8)
#' distribution(m,"X") <- binomial.lvm()
#' regression(m,latent.lifetime~X) <- 0.7
#' set.seed(30)
#' d <- sim(m,100)
#' fit.weib <- idm(formula02=Hist(observed.lifetime,event=seen.exit)~1,
#' formula01=Hist(time=list(L,R),event=seen.ill)~1,data=d,conf.int=FALSE)
#' plot(fit.weib)
#'
#' \dontrun{
#' ## FIXME: the limits for the 01 transition are a bit wide!?
#' ## with bootstrap confidence limits
#' fit.weib <- idm(formula02=Hist(observed.lifetime,event=seen.exit)~1,
#' formula01=Hist(time=list(L,R),event=seen.ill)~1,data=d,conf.int=TRUE)
#' plot(fit.weib)
#' }
#'  
#'
#' @S3method plot idm
plot.idm <- function(x,
                     conf.int=FALSE,
                     citype="shadow",
                     add=FALSE,
                     axes=TRUE,
                     col,
                     lwd,
                     lty,
                     xlim,
                     ylim,
                     xlab,
                     ylab,
                     legend=TRUE,
                     transition=c("01","02","12"),
                     ...){ 

    # {{{ collecting the (X, Y)-values of the lines
    if ((NCOL(x$time))>1){
        X01 <- x$time[,1]
        X02 <- x$time[,2]
        X12 <- x$time[,3]
        X <- unique(c(X01,X02,X12))
    }
    else {
        X01 <- X02 <- X12 <- X <- x$time
    }
    if(sum(c("01","02","12") %in% transition)!=length(transition))
        stop("'transition' must be a subset of 'c('01','02','12')'")
    Y <- list("01"=x$intensity01,"02"=x$intensity02,"12"=x$intensity12)
    Y <- Y[transition]
    nlines <- length(Y)
    # }}}
    # {{{ setting default arguments for plot, axes, legend, confint 
    if (missing(ylab)) ylab <- "Transition intensity"
    if (missing(xlab)) xlab <- "Time"
    if (missing(xlim)) xlim <- c(0, max(X))
    if (missing(ylim)) ylim <- c(0, 1)
    if (missing(lwd)) lwd <- rep(3,nlines)
    if (missing(col)) col <- 1:nlines
    if (missing(lty)) lty <- rep(1, nlines)
    if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
    if (length(lty) < nlines) lty <- rep(lty, nlines)
    if (length(col) < nlines) col <- rep(col, nlines)
    axis1.DefaultArgs <- list()
    axis2.DefaultArgs <- list(at=seq(0,1,.25),side=2,percent=TRUE)
    lines.DefaultArgs <- list(type="l")
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = xlab,ylab = ylab)
    legend.DefaultArgs <- list(legend=paste("Transition",names(Y)),lwd=lwd,col=col,lty=lty,cex=1.5,bty="n",y.intersp=1.3,x="topleft")
    confint.DefaultArgs <- list(x=x,citype="shadow",times=X,density=55,col=col[1:nlines],lwd=rep(2,nlines),lty=rep(3,nlines))
    # }}}
    control <- prodlim::SmartControl(call=  list(...),
                            keys=c("plot","lines","legend","confint","axis1","axis2"),
                            ignore=c("x","transition","add","col","lty","lwd","ylim","xlim","xlab","ylab","legend","conf.int","axes"),
                            defaults=list("plot"=plot.DefaultArgs,"lines"=lines.DefaultArgs,"legend"=legend.DefaultArgs,"confint"=confint.DefaultArgs,"axis1"=axis1.DefaultArgs,"axis2"=axis2.DefaultArgs),
                            forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),
                            ignore.case=TRUE,
                            replaceDefaults=FALSE,
                            verbose=TRUE)
    # {{{  plot and backGround
    if (!add) {
        do.call("plot",control$plot)
    }
    # }}}
    # {{{  axes
    if (!add) {
        if (axes){
            do.call("axis",control$axis1)
            if (control$axis2$percent & is.null(control$axis2$labels))
                control$axis2$labels <- paste(100*control$axis2$at,"%")
            do.call("axis",control$axis2[-match("percent",names(control$axis2),nomatch=0)])
        }
    }
    # }}}
    # {{{confidence intervals
    nix <- lapply(1:nlines,function(i){
        ci.lower <- x[[paste("lowerIntensity",names(Y)[i],sep="")]]
        ci.upper <- x[[paste("upperIntensity",names(Y)[i],sep="")]]
        time <- switch(i, "1"= X01,"2"=X02,"3"=X12)
        switch(citype,
               "bars"={
                   segments(x0=time,x1=time,y0=ci.lower,y1=ci.upper,lwd=lwd[i],col=col[i],lty=lty[i],...)
               },
               "shadow"={
                   ccrgb=as.list(col2rgb(col[i],alpha=TRUE))
                   names(ccrgb) <- c("red","green","blue","alpha")
                   ccrgb$alpha=control$confint$density
                   cc=do.call("rgb",c(ccrgb,list(max=255)))
                   ttt <- time
                   nt <- length(ttt)
                   ttt <- c(ttt,ttt)
                   uuu <- c(0,ci.upper[-nt],ci.upper)
                   lll <- c(0,ci.lower[-nt],ci.lower)
                   neworder <- order(ttt)
                   uuu <- uuu[neworder]
                   lll <- lll[neworder]
                   ttt <- sort(ttt)
                   polygon(x=c(ttt,rev(ttt)),y=c(lll,rev(uuu)),col=cc,border=NA)
               },{
                   lines(x=time,ci.lower,type="l",lwd=lwd[i],col=col[i],lty=3,...)
                   lines(x=time,ci.upper,type="l",lwd=lwd[i],col=col[i],lty=3,...)
               })
    })
    # }}}
    # {{{  adding the lines
    lines.type <- control$lines$type
    nix <- lapply(1:nlines, function(s) {
        time <- switch(s, "1"= X01,"2"=X02,"3"=X12)
        lines(x = time,y = Y[[s]],type = lines.type,col = col[s],lty = lty[s],lwd = lwd[s])
    })
    # {{{  legend
    if(legend==TRUE && !add && !is.null(names(Y))){
        if (is.null(control$legend$title)){
            if (x$method=="Splines")
                control$legend$title <- "M-spline intensity model"
            else
                control$legend$title <- "Weibull model"
        }
        save.xpd <- par()$xpd
        par(xpd=TRUE)
        do.call("legend",control$legend)
        par(xpd=save.xpd)
    }
    # }}}
    invisible(x)
}
# }}}
## if(class(x) == "idmWeib"){
## if (missing(transition)){
## if (missing(ylim)){
## ylim <- c(min(x$intensity01,x$intensity02,x$intensity12),max(x$intensity01,x$intensity02,x$intensity12))
## }
## if(conf.int){
## matplot(x$time, cbind(x$intensity01,x$lowerIntensity01,x$upperIntensity01),type="l",lty=c(1,1,1),col=c(1,1,1),lwd=c(2,1,1), xlab="Time",ylab="Weibull transition intensities",ylim=ylim)
## matlines(x$time, cbind(x$intensity02,x$lowerIntensity02,x$upperIntensity02),ylim=ylim, type="l",lwd=c(2,1,1),lty=c(2,2,2),col=c(2,2,2))
## matlines(x$time, cbind(x$intensity12,x$lowerIntensity12,x$upperIntensity12),ylim=ylim, type="l",lwd=c(2,1,1),lty=c(3,3,3),col=c(3,3,3))
## }else{
## matplot(x$time, cbind(x$intensity01,x$intensity02,x$intensity12),type="l",lwd=2,lty=c(1,2,3),col=c(1,2,3), xlab="Time",ylab="Weibull transition intensities",ylim=ylim, main=main)
## }	
## legend(pos.legend,c("01","02","12"),lty=c(1,2,3),lwd=c(2,2,2),col=c(1,2,3))
## }else{
## if(length(transition)==3){
## par(mfrow=c(3,1))
## if("01" %in% transition){
## if (missing(ylim)){
## ylim <- c(min(c(x$intensity01,x$lowerIntensity01,x$upperIntensity01)),max(c(x$intensity01,x$lowerIntensity01,x$upperIntensity01)))
## }
## if(!conf.int){
## plot(x$time, x$intensity01,type="l",col=1,lwd=2, xlab="Time",ylab="Weibull transition intensity 01",main=main)
## }else{
## matplot(x$time, cbind(x$intensity01,x$lowerIntensity01,x$upperIntensity01), col=c(1,1,1),type="l",lty=c(1,1,1),lwd=c(2,1,1),ylim=ylim
## ,xlab="Time",ylab="Weibull transition intensity 01", main=main)	
## #							legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
## }
## }
## if("02"%in% transition){
## if (missing(ylim)){
## ylim <- c(min(c(x$intensity02,x$lowerIntensityd02,x$upperIntensity02)),max(c(x$intensity02,x$lowerIntensity02,x$upperIntensity02)))
## }
## if(!conf.int){
## plot(x$time, x$intensity02, type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Weibull transition intensity 02",main=main)
## }else{
## matplot(x$time, cbind(x$intensity02,x$lowerIntensity02,x$upperIntensity02),type="l",col=c(2,2,2),lwd=c(2,1,1),lty=c(2,2,2),ylim=ylim
## ,xlab="Time",ylab="Weibull transition intensity 02",main=main)	
## #							legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
## }
## }				
## if("12" %in% transition){
## if (missing(ylim)){
## ylim <- c(min(c(x$intensity12,x$lowerIntensity12,x$upperIntensity12)),max(c(x$intensity12,x$lowerIntensity12,x$upperIntensity12)))
## }
## if(!conf.int){
## plot(x$time, x$intensity12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Weibull transition intensity 12", main=main)
## }else{
## matplot(x$time, cbind(x$intensity12,x$lowerIntensity12,x$upperIntensity12),type="l",col=c(3,3,3),lwd=c(2,1,1),lty=c(3,3,3),ylim=ylim
## ,xlab="Time",ylab="Weibull transition intensity 12",main=main)	
## #							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
## }
## }
## }else{	
## if(length(transition)==2){
## par(mfrow=c(2,1))
## if("01" %in% transition){
## if (missing(ylim)){
## ylim <- c(min(x$intensity01,x$lowerIntensity01,x$upperIntensity01),max(x$intensity01,x$lowerIntensity01,x$upperIntensity01))
## }
## if(!conf.int){
## plot(x$time, x$intensity01,type="l",col=1,lwd=2, xlab="Time",ylab="Weibull transition intensity 01",main=main)
## }else{
## matplot(x$time, cbind(x$intensity01,x$lowerIntensity01,x$upperIntensity01), col=c(1,1,1),type="l",lty=c(1,1,1),lwd=c(2,1,1),ylim=ylim
## ,xlab="Time",ylab="Weibull transition intensity 01", main=main)	
## #							legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
## }
## }
## if("02"%in% transition){
## if (missing(ylim)){
## ylim <- c(min(c(x$intensity02,x$lowerIntensity02,x$upperIntensity02)),max(c(x$intensity02,x$lowerIntensity02,x$upperIntensity02)))
## }
## if(!conf.int){
## plot(x$time, x$intensity02, type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Weibull transition intensity 02",main=main)
## }else{
## matplot(x$time, cbind(x$intensity02,x$lowerIntensity02,x$upperIntensity02),type="l",col=c(2,2,2),lwd=c(2,1,1),lty=c(2,2,2),ylim=ylim
## ,xlab="Time",ylab="Weibull transition intensity 02",main=main)	
## #							legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
## }
## }				
## if("12" %in% transition){
## if (missing(ylim)){
## ylim <- c(min(c(x$intensity12,x$lowerIntensity12,x$upperIntensity12)),max(c(x$intensity12,x$lowerIntensity12,x$upperIntensity12)))
## }
## if(!conf.int){
## plot(x$time, x$intensity12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Weibull transition intensity 12", main=main)
## }else{
## matplot(x$time, cbind(x$intensity12,x$lowerIntensity12,x$upperIntensity12),type="l",col=c(3,3,3),lwd=c(2,1,1),lty=c(3,3,3),ylim=ylim
## ,xlab="Time",ylab="Weibull transition intensity 12",main=main)	
## #							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
## }
## }
## }else{
## if(transition=="01"){
## if (missing(ylim)){
## ylim <- c(min(x$intensity01,x$lowerIntensity01,x$upperIntensity01),max(x$intensity01,x$lowerIntensity01,x$upperIntensity01))
## }
## if(!conf.int){
## plot(x$time,x$intensity01,type="l",col=1,lwd=2,xlab="Time",ylab="Weibull transition intensity 01",main=main)
## }else{
## matplot(x$time, cbind(x$intensity01,x$lowerIntensity01,x$upperIntensity01),type="l",col=c(1,1,1),lty=c(1,1,1),lwd=c(2,1,1),ylim=ylim
## ,xlab="Time",ylab="Weibull transition intensity 01", main=main)
## #							 legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,1,1),cex=0.7)
## }
## }
## if(transition=="02"){
## if (missing(ylim)){
## ylim <- c(min(x$intensity02,x$lowerIntensity02,x$upperIntensity02),max(x$intensity02,x$lowerIntensity02,x$upperIntensity02))
## }
## if(!conf.int){
## plot(x$time, x$intensity02,type="l",col=2,lwd=2,lty=2,xlab="Time",ylab="Weibull transition intensity 02",main=main)
## }else{
## matplot(x$time,cbind(x$intensity02,x$lowerIntensity02,x$upperIntensity02),type="l",col=c(2,2,2),lwd=c(2,1,1),lty=c(2,2,2),ylim=ylim
## ,xlab="Time",ylab="Weibull transition intensity 02",main=main)	
## #							 legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2),cex=0.7)	
## }
## }				
## if(transition=="12"){
## if (missing(ylim)){
## ylim <- c(min(x$intensity12,x$lowerIntensity12,x$upperIntensity12),max(x$intensity12,x$lowerIntensity12,x$upperIntensity12))
## }
## if(!conf.int){
## plot(x$time, x$intensity12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Weibull transition intensity 12", main=main)
## }else{
## matplot(x$time, cbind(x$intensity12,x$lowerIntensity12,x$upperIntensity12),type="l",col=c(3,3,3),lwd=c(2,1,1),lty=c(3,3,3),ylim=ylim
## ,xlab="Time",ylab="Weibull transition intensity 12", main=main)	
## #							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3),cex=0.7)
## }
## }
## }
## }	
## }
## }else{
## if (missing(transition)){
## xmin <- min(x$time)
## xmax <- max(x$time)
## if (missing(ylim)){
## ymin <- min(x$intensity01,x$intensity02,x$intensity12)
## ymax <- max(x$intensity01,x$intensity02,x$intensity12)
## ymin01 <- min(x$intensity01,x$lowerIntensity01,x$upperIntensity01)
## ymax01 <- max(x$intensity01,x$lowerIntensity01,x$upperIntensity01)
## ymin02 <- min(x$intensity02,x$lowerIntensity02,x$lowerIntensity02)
## ymax02 <- max(x$intensity02,x$lowerIntensity02,x$lowerIntensity02)
## ymin12 <- min(x$intensity12,x$lowerIntensity12,x$lowerIntensity12)
## ymax12 <- max(x$intensity12,x$lowerIntensity12,x$lowerIntensity12)
## ylim <- c(ymin,ymax)
## }
## if(conf.int){
## matplot(x$time[,1], cbind(x$intensity01,x$lowerIntensity01,x$upperIntensity01),type="l",lwd=c(2,1,1),lty=c(1,1,1),col=c(1,1,1),
## xlab="Time",ylab="Splines transition intensities",xlim=c(xmin,xmax),ylim=ylim, main=main)
## matlines(x$time[,2], cbind(x$intensity02,x$lowerIntensity02,x$upperIntensity02),xlim=c(xmin,xmax),ylim=ylim,type="l",lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2))
## matlines(x$time[,3], cbind(x$intensity12,x$lowerIntensity12,x$upperIntensity12),xlim=c(xmin,xmax),ylim=ylim,type="l",lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3))
## }else{
## plot(x$time[,1], x$intensity01, col=1, type="l",lwd=2, xlab="Time",ylab="Splines transition intensities",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
## main=main,)
## lines(x$time[,2], x$intensity02, col=2, type="l",lty=2,lwd=2,xlim=c(xmin,xmax),ylim=ylim)
## lines(x$time[,3], x$intensity12, col=3, type="l",lty=3,lwd=2,xlim=c(xmin,xmax),ylim=ylim)
## }
## legend(pos.legend,c("01","02","12"),lty=c(1,2,3),lwd=c(2,2,2),col=c(1,2,3))	
## }else{
## if(length(transition)==3){
## par(mfrow=c(3,1))
## if("01" %in% transition){
## if (missing(ylim)){
## ylim <- c(min(x$intensity01,x$lowerIntensity01,x$upperIntensity01),max(x$intensity01,x$lowerIntensity01,x$upperIntensity01))
## }
## if(!conf.int){
## plot(x$time[,1],x$intensity01,type="l",col=1,lwd=2, xlab="Time",ylab="Splines transition intensity 01", main=main)
## }else{
## matplot(x$time[,1],cbind(x$intensity01,x$lowerIntensity01,x$upperIntensity01),type="l",col=c(1,1,1),lwd=c(2,1,1),lty=c(1,1,1),ylim=ylim
## ,xlab="Time",ylab="Splines transition intensity", main=main)
## #							 legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,2,3),cex=0.7)
## }
## }
## if("02"%in% transition){
## if (missing(ylim)){
## ylim <- c(min(c(x$intensity02,x$lowerIntensity02,x$upperIntensity02)),max(c(x$intensity02,x$lowerIntensity02,x$upperIntensity02)))
## }
## if(!conf.int){
## plot(x$time[,2],x$intensity02,type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Splines transition intensity 02", main=main)
## }else{
## matplot(x$time[,2],cbind(x$intensity02,x$lowerIntensity02,x$upperIntensity02),type="l",col=c(2,2,2),lty=c(2,2,2),lwd=c(2,1,1),ylim=ylim
## ,xlab="Time",ylab="Splines transition intensity", main=main)
## #							 legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2),cex=0.7)
## }
## }				
## if("12" %in% transition){
## if (missing(ylim)){
## ylim <- c(min(c(x$intensity12,x$lowerIntensity12,x$upperIntensity12)),max(c(x$intensity12,x$lowerIntensity12,x$upperIntensity12)))
## }
## if(!conf.int){
## plot(x$time[,3],x$intensity12,type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Splines transition intensity 12", main=main)
## }else{
## matplot(x$time[,3],cbind(x$intensity12,x$lowerIntensity12,x$upperIntensity12),type="l",col=c(3,3,3),lty=c(3,3,3),lwd=c(2,1,1),ylim=ylim
## ,xlab="Time",ylab="Splines transition intensity", main=main)
## #							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3),cex=0.7)
## }
## }
## }else{
## if(length(transition)==2){
## par(mfrow=c(2,1))
## if("01" %in% transition){
## if (missing(ylim)){
## ylim <- c(min(x$intensity01,x$lowerIntensity01,x$upperIntensity01),max(x$intensity01,x$lowerIntensity01,x$upperIntensity01))
## }
## if(!conf.int){
## plot(x$time[,1],x$intensity01, type="l",col=1,lwd=2, xlab="Time",ylab="Splines transition intensity 01", main=main)
## }else{
## matplot(x$time[,1],cbind(x$intensity01,x$lowerIntensity01,x$upperIntensity01),type="l",col=c(1,1,1),lwd=c(2,1,1),lty=c(1,1,1),ylim=ylim
## ,xlab="Time",ylab="Splines transition intensity", main=main)
## #							 legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,1,1),cex=0.7)
## }
## }
## if("02"%in% transition){
## if (missing(ylim)){
## ylim <- c(min(x$intensity02,x$lowerIntensity02,x$upperIntensity02),max(x$intensity02,x$lowerIntensity02,x$upperIntensity02))
## }
## if(!conf.int){
## plot(x$time[,2],x$intensity02, type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Splines transition intensity 02", main=main)
## }else{
## matplot(x$time[,2],cbind(x$intensity02,x$lowerIntensity02,x$upperIntensity02),type="l",col=c(2,2,2),lty=c(2,2,2),lwd=c(2,1,1),ylim=ylim
## ,xlab="Time",ylab="Splines transition intensity", main=main)
## #							 legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2),cex=0.7)
## }
## }				
## if("12" %in% transition){
## if (missing(ylim)){
## ylim <- c(min(x$intensity12,x$lowerIntensity12,x$upperIntensity12),max(x$intensity12,x$lowerIntensity12,x$upperIntensity12))
## }
## if(!conf.int){
## plot(x$time[,3],x$intensity12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Splines transition intensity 12", main=main)
## }else{
## matplot(x$time[,3],cbind(x$intensity12,x$lowerIntensity12,x$upperIntensity12),type="l",col=c(3,3,3),lty=c(3,3,3),lwd=c(2,1,1),ylim=ylim
## ,xlab="Time",ylab="Splines transition intensity", main=main)
## #							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3),cex=0.7)
## }
## }					
## }else{
## if(transition=="01"){
## if (missing(ylim)){
## ylim <- c(min(x$intensity01,x$lowerIntensity01,x$upperIntensity01),max(x$intensity01,x$lowerIntensity01,x$upperIntensity01))
## }
## if(!conf.int){
## plot(x$time[,1],x$intensity01,type="l",col=1,lwd=2,xlab="Time",ylab="Splines transition intensity 01",main=main)
## }else{
## matplot(x$time[,1],cbind(x$intensity01,x$lowerIntensity01,x$upperIntensity01),type="l",col=c(1,1,1),lwd=c(2,1,1),lty=c(1,1,1),ylim=ylim
## ,xlab="Time",ylab="Splines transition intensity", main=main)
## #							 legend(pos.legend,c("01","lower","upper"),lty=c(1,1,1),lwd=c(2,1,1),col=c(1,1,1),cex=0.7)
## }
## }
## if(transition=="02"){
## if (missing(ylim)){
## ylim <- c(min(x$intensity02,x$lowerIntensity02,x$upperIntensity02),max(x$intensity02,x$lowerIntensity02,x$upperIntensity02))
## }
## if(!conf.int){
## plot(x$time[,2],x$intensity02, type="l",col=2,lwd=2,lty=2, xlab="Time",ylab="Splines transition intensity 02", main=main)
## }else{
## matplot(x$time[,2],cbind(x$intensity02,x$lowerIntensity02,x$upperIntensity02),type="l",col=c(2,2,2),lty=c(2,2,2),lwd=c(2,1,1),ylim=ylim
## ,xlab="Time",ylab="Splines transition intensity", main=main)
## #							 legend(pos.legend,c("02","lower","upper"),lty=c(2,2,2),lwd=c(2,1,1),col=c(2,2,2),cex=0.7)
## }
## }				
## if(transition=="12"){
## if (missing(ylim)){
## ylim <- c(min(x$intensity02,x$lowerIntensity02,x$upperIntensity02),max(x$intensity02,x$lowerIntensity02,x$upperIntensity02))
## }
## if(!conf.int){
## plot(x$time[,3],x$intensity12, type="l",col=3,lwd=2,lty=3, xlab="Time",ylab="Splines transition intensity 12", main=main)
## }else{
## matplot(x$time[,3],cbind(x$intensity12,x$lowerIntensity12,x$upperIntensity12),type="l",col=c(3,3,3),lty=c(3,3,3),lwd=c(2,1,1),ylim=ylim
## ,xlab="Time",ylab="Splines transition intensity", main=main)
## #							 legend(pos.legend,c("12","lower","upper"),lty=c(3,3,3),lwd=c(2,1,1),col=c(3,3,3),cex=0.7)
## }
## }
## }
## }
## }
## }
## return(invisible())
## }
