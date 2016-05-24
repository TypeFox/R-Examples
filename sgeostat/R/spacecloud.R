#
# spacecloud.s creates a box plot of squared or square root difference plots...
assign("spacecloud",
function(point.obj,pair.obj,a1,a2,type='r',query.a=NULL,...) {

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

  plot(pair.obj$dist,diff,
          xlab='lag',ylab=ylab,...)
  title(deparse(substitute(point.obj)))
  if(!is.null(query.a)) {
    query.att <- point.obj[[match(query.a,names(point.obj))]]
    cat('Identify "from" points...')
    identify(pair.obj$dist,diff,
             query.att[pair.obj$from])#,col=2)
    cat('\nIdentify "to" points...')
    identify(pair.obj$dist,diff,
             query.att[pair.obj$to])#,col=3)

    cat('\n')
  }

  return(invisible(NULL))
})
