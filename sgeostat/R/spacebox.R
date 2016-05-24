# spacebox.s creates a box plot of squared or square root difference plots...
assign("spacebox",
function(point.obj,pair.obj,a1,a2,type='r') {

  if (!inherits(point.obj,"point")) stop('Point.obj must be of class, "point".\n')

  if (!inherits(pair.obj,"pair")) stop('Pair.obj must be of class, "pair".\n')

  if(missing(a1)) stop('Must enter at least one attribute.\n')
  if(missing(a2)) a2 <- a1

  a1 <- point.obj[[match(a1,names(point.obj))]]

  a2 <- point.obj[[match(a2,names(point.obj))]]

  if (type=='r') {  # square root difference cloud
    diff <- (abs(a1[pair.obj$from]-a2[pair.obj$to]))^0.5
    ylab <- 'square root differnece'
  }
  else {
    diff <- (a1[pair.obj$from]-a2[pair.obj$to])^2
    ylab <- 'squared differnece'
  }
  names<-sort(unique(pair.obj$lags))
  names <- levels(pair.obj$lags)
  boxplot(split(diff,pair.obj$lags),
# revision 11/4-99 rsb
#	  names=names,
	  names.x=names,
	  varwidth=TRUE)
  title(xlab='lag',
	ylab=ylab)

  return(invisible(NULL))
})
