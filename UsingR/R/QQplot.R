##' make a fancy qqplot
##'
##' @param x a data vector
##' @param y a data vector
##' @param n split into this many quantile regions
##' @param xsf scale factor
##' @param ysf scale factor
##' @param main title
##' @param xlab x label
##' @param ylab y label
##' @param pch plot character
##' @param pcol point color
##' @param shade shade color
##' @param ... passed along
##'
##' @return NULL
##' @export
QQplot = function(x, y, n = 20,
  xsf= 4,ysf= 4,
  main="qqplot",xlab=deparse(substitute(x)), ylab = deparse(substitute(y)),
  pch =16,pcol = "black",                # points
  shade = "gray",
  ...) {
  ## n -- number of points
  ## xsf scale factor for x density
  ## ysf same for y

  ## scale x and y
  x.sc = (x-mean(x))/sd(x); y.sc = (y-mean(y))/sd(y)
  
  zip = function(x,y,do.na = TRUE) {
    ## alternate x,y
    k = length(x)
    ret = c()
    m = 2
    if(do.na == TRUE)
      m = 3
    ret[m*(0:(k-1))+1] = x
    ret[m*(0:(k-1))+2] = y
    if(m == 3)
      ret[m*(0:(k-1))+3] = rep(NA,k)             # recycle
    return(ret)
  }
  zip4 = function(x,y,z,w) {
    ## alternate x,y,z,w
    k = length(x)
    ret = c()
    m = 4
    ret[m*(0:(k-1))+1] = x
    ret[m*(0:(k-1))+2] = y
    ret[m*(0:(k-1))+3] = z
    ret[m*(1:k)] = w
#    ret[m*(1:k)] = rep(NA,k)
    return(ret)
  }

  ## return y value for given x value
  findInDensity= function(x,dens) {
    if(class(dens) == "density" & !is.na(x)) {
      if(any(dens$x < x)) {
        i = max(which(dens$x < x))
        dens$y[i]
      } else {
        dens$y[length(dens$y)]
      }
    } else {
      return(NA)
    }
  }
  
  
  ## n should be even
  k = n %/% 2                           # half of n
  n = k * 2 - 1
  vals = (1:(n-1))/(n)
  
  dens = list(
    x= density(x.sc),
    y = density(y.sc)
    )
  qpts = list(
    x = quantile(x.sc,vals),
    y = quantile(y.sc,vals)
    )


  xlim = c(min(dens$x$x) - ysf*max(dens$y$y),max(dens$x$x))
  ylim = c(min(dens$y$x) - xsf*max(dens$x$y),max(dens$y$x))

  
  plot.new()
  plot.window(xlim=xlim,ylim=ylim,...)
  title(main=main,xlab=xlab,ylab=ylab)
  
  lines(xlim,rep(min(dens$y$x),2))
  lines(rep(min(dens$x$x),2),ylim)
  
  ## x first
  lines(dens$x$x,min(dens$y$x) - xsf*dens$x$y, type="l",xaxt="n",yaxt="n")

  ## shade in.
  new.x = zip4(qpts$x[2*(1:k)-1],qpts$x[2*(1:k)-1],
    qpts$x[2*(1:k)],qpts$x[2*(1:k)])
  new.y = zip4(rep(0,k),
    -sapply(qpts$x[2*(1:k)-1],function(x) findInDensity(x,dens$x)),
    -sapply(qpts$x[2*(1:k)],function(x) findInDensity(x,dens$x)),
    rep(0,k))
  polygon(c(min(dens$x$x),new.x,max(dens$x$x)),min(dens$y$x)+xsf*c(0,new.y,0), col=shade)

  retval = list(
    x=new.x,
    y=new.y)

  
  ## y now
  lines(min(dens$x$x)-ysf*dens$y$y,dens$y$x,xaxt="n",yaxt="n")

  ## shade in.
  new.x = zip4(rep(0,k),
    -sapply(qpts$y[2*(1:k)-1],function(x) findInDensity(x,dens$y)),
    -sapply(qpts$y[2*(1:k)],  function(x) findInDensity(x,dens$y)),
    rep(0,k))
  new.y = zip4(qpts$y[2*(1:k)-1],qpts$y[2*(1:k)-1],
    qpts$y[2*(1:k)],qpts$y[2*(1:k)])
  polygon(min(dens$x$x)+ysf*c(new.x,max(new.x)),c(new.y,0), col=shade)

  ## plot densities

  points(qpts$x,qpts$y, pch=pch, col=pcol)


  lty = 2
  lines(zip(qpts$x,qpts$x),zip(rep(min(dens$y$x),n-1),qpts$y), lty=lty)
  lines(zip(rep(min(dens$x$x),n-1),qpts$x),zip(qpts$y,qpts$y), lty=lty)
}


