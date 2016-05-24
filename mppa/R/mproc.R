##setClass("mproc", representation(events = "list", start="numeric", end="numeric", mks="list", type="list"))
setClass("mproc", representation(events = "list", start="numeric", end="numeric"))


setMethod("plot", signature(x="mproc"), function(x, period=NA,id=1,mks=NULL,palette=rainbow, cols=NULL,xlab="Time", ylab="",...){
  start=ifelse(length(x@start)==0, min(unlist(lapply(x@events, function(x){min(as.numeric(x))}))), x@start)
  end=ifelse(length(x@end)==0, max(unlist(lapply(x@events, function(x){max(as.numeric(x))}))), x@end)

  if (is.na(period)){
    x = x@events
    plot.default(x=c(), y=c(), yaxt="n", xlim=c(start, end), xlab=xlab, ylab=ylab, ylim = c(0, length(x)+1),...)
    if (is.null(names(x)))ns=1:length(x)
    else ns = names(x)
    axis(2, at=length(x):1, labels=ns)
    for (i in 1:length(x)){
      if (!is.null(cols)){
        if (!is.list(cols)){cols=list(cols)}
        points(as.numeric(x[[i]]), rep(length(x)-i+1, length(x[[i]])), col=cols[[i]],...)
      }
      else if (!is.null(mks)){
        if (!is.list(mks)){mks=list(mks)}
        ##browser()
        colindex=tapply(x[[i]],mks[[i]]); ncols=max(colindex)
        points(as.numeric(x[[i]]), rep(length(x)-i+1, length(x[[i]])), col=palette(ncols)[colindex],...)
      }
      else {
        points(as.numeric(x[[i]]), rep(length(x)-i+1, length(x[[i]])), ...)
      }
    }    
  }
  
  else{
    if (ylab=="")ylab="Days"
    #browser()
    x = as.numeric(x@events[[id]])
    d=sapply(x, function(e){
      day=(e-start)%/%period+1; time = (e-start)%%period
      c(time,day)
    })
    fulldays = (end-start)%/%period
    if (!is.null(cols)){
      if (!is.list(cols)){cols=list(cols)}
      plot(d[1,], (fulldays+1)-d[2,]+1, xlim=c(0,period), ylim=c(1,fulldays+1), yaxt="n", xlab=xlab, ylab=ylab,type='p',col=cols[[id]],...)     
    }
    else if (!is.null(mks)){
      if (!is.list(mks)){mks=list(mks)}
      if (length(mks)==1){mks=mks[[1]]} else mks=mks[[id]]
      colindex=tapply(x,mks); ncols=max(colindex)
      plot(d[1,], (fulldays+1)-d[2,]+1, xlim=c(0,period), ylim=c(1,fulldays+1), yaxt="n", xlab=xlab, ylab=ylab,type='p',col=palette(ncols)[colindex],...)
    }
    else{
      plot(d[1,], (fulldays+1)-d[2,]+1, xlim=c(0,period), ylim=c(1,fulldays+1), yaxt="n", xlab=xlab, ylab=ylab,type='p',...)
    }
    #browser()
    axis(2, at=1:(fulldays+1), labels=(fulldays+1):1)
  }
})

## setMethod("plot", signature(x="mproc", y="missing"), function(x, period=NA,id=1,mks=NULL,palette=rainbow, cols=NULL,ylab="", xlab="Time",...){
##   start=ifelse(length(x@start)==0, min(unlist(lapply(x@events, function(x){min(as.numeric(x))}))), x@start)
##   end=ifelse(length(x@end)==0, max(unlist(lapply(x@events, function(x){max(as.numeric(x))}))), x@end)

##   if (is.na(period)){
##     x = x@events
##     plot.default(x=c(), y=c(), yaxt="n", xlim=c(start, end), xlab=xlab, ylab=ylab, ylim = c(0, length(x)+1),...)
##     if (is.null(names(x)))ns=1:length(x)
##     else ns = names(x)
##     axis(2, at=length(x):1, labels=ns)
##     for (i in 1:length(x)){
##       if (!is.null(cols)){
##         if (!is.list(cols)){cols=list(cols)}
##         points(as.numeric(x[[i]]), rep(length(x)-i+1, length(x[[i]])), col=cols[[i]],...)
##       }
##       else if (!is.null(mks)){
##         if (!is.list(mks)){mks=list(mks)}
##         ##browser()
##         colindex=tapply(x[[i]],mks[[i]]); ncols=max(colindex)
##         points(as.numeric(x[[i]]), rep(length(x)-i+1, length(x[[i]])), col=palette(ncols)[colindex],...)
##       }
##       else {
##         points(as.numeric(x[[i]]), rep(length(x)-i+1, length(x[[i]])), ...)
##       }
##     }    
##   }
  
##   else{
##     if (ylab=="")ylab="Days"
##     #browser()
##     x = as.numeric(x@events[[id]])
##     d=sapply(x, function(e){
##       day=(e-start)%/%period+1; time = (e-start)%%period
##       c(time,day)
##     })
##     fulldays = (end-start)%/%period
##     if (!is.null(cols)){
##       if (!is.list(cols)){cols=list(cols)}
##       plot(d[1,], (fulldays+1)-d[2,]+1, xlim=c(0,period), ylim=c(1,fulldays+1), yaxt="n", xlab=xlab, ylab=ylab,type='p',col=cols[[id]],...)     
##     }
##     else if (!is.null(mks)){
##       if (!is.list(mks)){mks=list(mks)}
##       if (length(mks)==1){mks=mks[[1]]} else mks=mks[[id]]
##       colindex=tapply(x,mks); ncols=max(colindex)
##       plot(d[1,], (fulldays+1)-d[2,]+1, xlim=c(0,period), ylim=c(1,fulldays+1), yaxt="n", xlab=xlab, ylab=ylab,type='p',col=palette(ncols)[colindex],...)
##     }
##     else{
##       plot(d[1,], (fulldays+1)-d[2,]+1, xlim=c(0,period), ylim=c(1,fulldays+1), yaxt="n", xlab=xlab, ylab=ylab,type='p',...)
##     }
##     #browser()
##     axis(2, at=1:(fulldays+1), labels=(fulldays+1):1)
##   }
## })

place <- function(e, m, period, ...){
  start=m@start; end=m@end
  day=(e-start)%/%period+1; time = (e-start)%%period
  fulldays = (end-start)%/%period
  points(time, (fulldays+1)-day+1, ...)
}


mp <- function(..., start=numeric(0), end=numeric(0)){  
  if (length(list(...))>1){events=list(...)}
  else if (is.list(...)){events=(...)}
  else {events=list(...)}    
  start=ifelse(length(start)==0, min(unlist(lapply(events, function(x){min(as.numeric(x))}))), start)
  end=ifelse(length(end)==0, max(unlist(lapply(events, function(x){max(as.numeric(x))}))),end)
  new("mproc", events=events, start=start, end=end)
}

## mp <- function(..., start=numeric(0), end=numeric(0), mks=list(), type=list()){
##   if (length(list(...))>1){events=list(...)}
##   else if (is.list(...)){events=(...)}
##   else {events=list(...)}    
##   start=ifelse(length(start)==0, min(unlist(lapply(events, function(x){min(as.numeric(x))}))), start)
##   end=ifelse(length(end)==0, max(unlist(lapply(events, function(x){max(as.numeric(x))}))),end)
##   new("mproc", events=events, start=start, end=end, mks=mks, type=type)
## }
