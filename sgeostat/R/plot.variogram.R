################
#  Create a plotting routine for variograms...
assign("plot.variogram",	
function(x,var.mod.obj=NULL,title.str=NULL,ylim,type='c',N=FALSE,...) {

#  oldpar <- par()
#  par(mfrow=c(2,1),lab=c(12,5,7), #lab=c(length(x$lag),5,7),
#      mar=c(4,12,4,12)+.1)

if (!inherits(x,'variogram')) stop('x must be of class "variogram"')
if (!missing(var.mod.obj))
  if (!inherits(var.mod.obj,'variogram.model')) stop('x must be of class "variogram.model"')

if(type!='c'&type!='r'&type!='m') stop('type must be "c", "r", or "m".\n')

  if(type=='c'){
    ylabel <- 'Classical semi-variogram estimator'
    y <- x$classic
  }
  if(type=='r'){
    ylabel <- 'Robust semi-variogram estimator'
    y <- x$robust
  }
  if(type=='m'){
    ylabel <- 'Median semi-varigram estimator'
    y <- x$med
  }
  y <- y/2


  if(missing(ylim)) ylim <- c(0,max(x$classic/2,x$robust/2,na.rm=TRUE))
  plot(x$bins,y,
       ylim=ylim,
       xlim=c(0,max(x$bins)),
       xlab="Lag",ylab=ylabel,
       type="p")

  if(N)
    text(x$bins,y,x$n)

  if(is.null(title.str))
    title(paste("Variogram estimator:",deparse(substitute(x))))
  else
    title(title.str)
				    
# See if we need to plot a fitted variogram...
  if(!is.null(var.mod.obj)) {
    if (is.null(attr(var.mod.obj,"class"))) stop('var.mod.obj must be of class, "variogram.model".\n')
    else if (attr(var.mod.obj,"class") != 'variogram.model') stop('var.mod.obj must be of class, "variogram.model".\n')
    h <- seq(from=0.0001,to=max(x$bins),length=50)
    lines(h,var.mod.obj$model(h,var.mod.obj$parameters))
  }

#  plot(x$lag,x$med,
#       ylim=c(0,max(x$classic,x$robust,x$med,na.rm=TRUE)),
#       xlab="Lag",ylab="Median estimator",
#       type="h")

#  plot(x$bins,x$robust,
#       ylim=c(0,max(x$classic,x$robust,x$med,na.rm=TRUE)),
#       xlim=c(0,max(x$bins)),
#       xlab="Lag",ylab="Robust estimator",
#       type="p")


#  par(mfrow=c(1,1))
#  invisible(par(oldpar))
})
