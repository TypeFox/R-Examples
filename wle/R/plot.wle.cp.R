#############################################################
#                                                           #
#	plot.wle.cp function                                    #
#	E-mail: claudio@unive.it                                #
#	Date: October, 27, 2003                                 #
#	Version: 0.4-2                                          #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli                  #
#                                                           #
#############################################################

plot.wle.cp <- function(x, base.line=0, num.max=20, plot.it=TRUE, log.scale=FALSE, xlab="Number of Predictors", ylab=NULL, verbose=FALSE, ...) {

    if (is.null(x$terms)) {
        stop("invalid \'wle.cp\' object")
    }

    wcp <- x$wcp

    if (num.max<1) {
        if (verbose) cat("plot.wle.cp: num.max can not less than 1, num.max=1 \n")
        num.max <- 1
    }

    if (is.null(nrow(wcp)) | nrow(wcp)==1) {
        num.model <- 1
    } else {
        num.model <- nrow(wcp) 
    }

    if (num.model<num.max) {
        if (verbose) cat("plot.wle.cp: The number of models is less than num.max \n")
        num.max <- num.model
    }

    if (is.null(ncol(wcp))) {
        stop("No models to plot")
    } else {
        nvar <- ncol(wcp)-1
    }

good.model <- (apply(wcp[,1:nvar],1,sum)+base.line>=wcp[,nvar+1])
wcp.good <- matrix(wcp[good.model,],ncol=nvar+1)
wcp.bad <- matrix(wcp[!good.model,],ncol=nvar+1)
ordine.good <- order(wcp.good[,nvar+1])
ordine.bad <- order(wcp.bad[,nvar+1])
wcp.good <- matrix(wcp.good[ordine.good,],ncol=(nvar+1))
num.good <- dim(wcp.good)[1]
wcp.bad <- matrix(wcp.bad[ordine.bad,],ncol=(nvar+1))
num.bad <- dim(wcp.bad)[1]

label.good <- character()
for(i in 1:nvar){
label.good <- paste(label.good,wcp.good[,i],sep="")
}
label.bad <- character()
for(i in 1:nvar){
label.bad <- paste(label.bad,wcp.bad[,i],sep="")
}


xcoord.good <- apply(matrix(wcp.good[,1:nvar],ncol=nvar),1,sum)[1:min(num.max,num.good)]
ycoord.good <- wcp.good[,nvar+1][1:min(num.max,num.good)]

label.good <- label.good[1:min(num.max,num.good)]

xcoord.best <- xcoord.good[1]
ycoord.best <- ycoord.good[1]

label.best <- label.good[1]

if (length(xcoord.good)==1) {
    xcoord.good <- xcoord.best <- 0
    ycoord.good <- ycoord.best <- 0
    plot.good <- FALSE
} else {
    xcoord.good <- xcoord.good[-1]
    ycoord.good <- ycoord.good[-1]
    label.good <- label.good[-1]
    plot.good <- TRUE
}

if (num.max>num.good) {
    xcoord.bad <- apply(matrix(wcp.bad[,1:nvar],ncol=nvar),1,sum)[1:min(num.bad,num.max-num.good)]
    ycoord.bad <- wcp.bad[,nvar+1][1:min(num.bad,num.max-num.good)]
    label.bad <- label.bad[1:min(num.bad,num.max-num.good)]
    plot.bad <- TRUE
} else {
    xcoord.bad <- 0
    ycoord.bad <- 0
    plot.bad <- FALSE
}

xlim.min <- min(xcoord.good, xcoord.bad, xcoord.best)
xlim.max <- max(xcoord.good, xcoord.bad, xcoord.best)

yetichetta <- "WCp"

if (log.scale) {
    ycoord.good <- log10(ycoord.good+min(ycoord.good,ycoord.bad,ycoord.best)+1)
    ycoord.bad <- log10(ycoord.bad+min(ycoord.good,ycoord.bad,ycoord.best)+1)
    ycoord.best <- log10(ycoord.best+min(ycoord.good,ycoord.bad,ycoord.best)+1)
    yetichetta <- "WCp log10 scale"
}

ylim.min <- min(ycoord.good, ycoord.bad, ycoord.best)
ylim.max <- max(ycoord.good, ycoord.bad, ycoord.best)

if (is.null(ylab)) {
    ylab <- yetichetta
}

if(plot.it)
{
plot(xcoord.best,ycoord.best,xlim=c(xlim.min,xlim.max),ylim=c(ylim.min,ylim.max),xlab=xlab,ylab=ylab,type="n")
text(xcoord.best,ycoord.best,col=4,labels=label.best)

if(plot.good)
{
text(xcoord.good,ycoord.good,col=3,labels=label.good)
}

if(plot.bad)
{
text(xcoord.bad,ycoord.bad,col=2,labels=label.bad)
}


if(!log.scale)
{
abline(base.line,1,col=2)
abline(0,1)
}
else
{
vettx <- seq(xlim.min,xlim.max,0.5)
vetty <- log10(vettx+min(ycoord.good,ycoord.bad,ycoord.best)+1)
vetty.base.line <- log10(vettx+min(ycoord.good,ycoord.bad,ycoord.best)+1+base.line)
lines(vettx,vetty.base.line,col=2,type="l")
lines(vettx,vetty,type="l")
}

}

invisible(list(num.good=num.good,num.bad=num.bad,wcp.good=wcp.good, wcp.bad=wcp.bad))
}
