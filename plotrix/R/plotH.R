plotH <- function (x,...) {
  UseMethod("plotH") 
}

plotH.formula <- function(x,data=NULL,xlab=names(mf)[2],ylab=names(mf)[1],...) {
  mf <- model.frame(x,data=data)                               # get model frame
  if (ncol(mf)>2) stop("Function currently only accepts one variable on RHS of formula")
  plotH.default(mf[,2],mf[,1],xlab=xlab,ylab=ylab,...)
}

plotH.default <- function(x,y,xlab=paste(deparse(substitute(x))),ylab=paste(deparse(substitute(y))),width=0.6,ylim=c(0,max(y)),col="gray",...) {
  plotHq <- function(x,y,xlab,ylab,width,ylim,col,...) {
    if (is.null(width)) width <- min(diff(x))
    xleft <- x-width/2
    xright <- x+width/2
    ytop <- y
    if (is.null(ylim)) {
      ybottom <- min(y[is.finite(y)])
      ylim <- c(ybottom,max(y[is.finite(y)]))
    } else ybottom <- min(ylim)
    plot(x,y,type="n",ylim=ylim,xlab=xlab,ylab=ylab,...)
    rect(xleft,ybottom,xright,ytop,col=col,...)
  }  # end plotHq internal function

  plotHc <- function(x,y,xlab,ylab,width,ylim,col,...) {
    names(y) <- x
    barplot(y,xlab=xlab,ylab=ylab,width=width,ylim=ylim,col=col,...)
  }  # end plotHc internal function
 ## Start of main function
  if (!is.numeric(y)) stop("Y (or LHS) variable must be quantitative.")
  if (is.numeric(x)) plotHq(x,y,xlab,ylab,width,ylim,col,...)
    else plotHc(x,y,xlab,ylab,width,ylim,col,...) 
}
