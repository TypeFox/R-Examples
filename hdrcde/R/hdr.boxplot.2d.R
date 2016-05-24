# Functions for 2-d HDR calculation
# Copyright Rob J Hyndman, 2013.


# Main function. Called by hdr.boxplot.2d

hdr.2d <- function(x, y, prob = c(50, 95, 99), den=NULL, kde.package=c("ash","ks"), h=NULL, 
    xextend=0.15, yextend=0.15)
{
   # Convert prob to coverage percentage if necessary
  if(max(prob) > 50) # Assume prob is coverage percentage
    alpha <- (100-prob)/100
  else # prob is tail probability (for backward compatibility)
    alpha <- prob
  alpha <- sort(alpha)

  # Estimate bivariate density
  if(is.null(den))
    den <- den.estimate.2d(x,y,kde.package,h,xextend,yextend)

  # Calculates falpha needed to compute HDR of bivariate density den.
  # Also finds approximate mode.
  fxy <- interp.2d(den$x,den$y,den$z,x,y)
  falpha <- quantile(fxy, alpha)
  index <- which.max(fxy)
  mode <- c(x[index],y[index])
  return(structure(list(mode=mode,falpha=falpha,fxy=fxy, den=den, alpha=alpha, x=x, y=y), class="hdr2d"))
}

hdr.boxplot.2d <- function(x, y, prob=c(50, 99), kde.package=c("ash","ks"), h=NULL,
  xextend=0.15, yextend=0.15, xlab="", ylab="", 
  shadecols=gray((length(prob):1)/(length(prob)+1)), pointcol=1, ...)
{
  # Estimate bivariate density
  hdr <- hdr.2d(x, y, prob=prob, kde.package=kde.package, h=h, xextend=xextend, yextend=yextend)
  # Produce plot
  plot(hdr, xlab=xlab, ylab=ylab, shadecols=shadecols, pointcol=pointcol, outside.points=TRUE,...)
}

plot.hdr2d <- function(x, shaded=TRUE, show.points=FALSE, outside.points=FALSE, pch=20,  shadecols=gray((length(x$alpha):1)/(length(x$alpha)+1)), 
    pointcol=1, ...)
{
  if(shaded)
    hdrcde.filled.contour(x$den$x,x$den$y,x$den$z,levels=c(x$falpha,1e10),col=shadecols,...)
  else
    contour(x$den,levels=x$falpha,labcex=0,...)
  if(show.points)
    points(x$x, x$y, pch=pch, col=pointcol)
  else if(outside.points)
  {
    index <- (x$fxy < 0.99999*min(x$falpha))
    points(x$x[index], x$y[index], pch=pch, col=pointcol)
  }
  points(x$mode[1],x$mode[2],pch=19,col="black")
  invisible(x)
}

hdrcde.filled.contour <- function (x,y,z, xlim = range(x, finite = TRUE), 
    ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
    levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
    col = color.palette(length(levels) - 1), plot.title, plot.axes, 
    asp = NA, xaxs = "i", yaxs = "i", las = 1, 
    axes = TRUE, frame.plot = axes, ...) 
{
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  .filled.contour(x,y,z,levels,col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
      box()
  if (missing(plot.title)) 
      title(...)
  else plot.title
  invisible()
}


"interp.2d" <- function(x, y, z, x0, y0)
{
  # Bilinear interpolation of function (x,y,z) onto (x0,y0).
  # Taken from Numerical Recipies (second edition) section 3.6.
  # Called by hdr.2d
  # Vectorized version of old.interp.2d. 
  # Contributed by Mrigesh Kshatriya (mkshatriya@zoology.up.ac.za)

  nx <- length(x)
  ny <- length(y)
  n0 <- length(x0)
  z0 <- numeric(length = n0)
  xr <- diff(range(x))
  yr <- diff(range(y))
  xmin <- min(x)
  ymin <- min(y)
  j <- ceiling(((nx - 1) * (x0 - xmin))/xr)
  k <- ceiling(((ny - 1) * (y0 - ymin))/yr)
  j[j == 0] <- 1
  k[k == 0] <- 1
  j[j == nx] <- nx - 1
  k[k == ny] <- ny - 1
  v <- (x0 - x[j])/(x[j + 1] - x[j])
  u <- (y0 - y[k])/(y[k + 1] - y[k]) 
  AA <- z[cbind(j, k)]
  BB <- z[cbind(j + 1, k)]
  CC <- z[cbind(j + 1, k + 1)]
  DD <- z[cbind(j, k + 1)]
  z0 <- (1 - v) * (1 - u) * AA + v * (1 - u) * BB + v * u * CC + (1 - v) * u * DD
  return(z0)
}


# Bivariate density estimate

den.estimate.2d <- function(x, y, kde.package=c("ash","ks"), h=NULL, xextend=0.15,yextend=0.15)
{
  kde.package <- match.arg(kde.package)
  # Find ranges for estimates
  xr <- diff(range(x,na.rm=TRUE))
  yr <- diff(range(y,na.rm=TRUE))
  xr <- c(min(x)-xr*xextend,max(x)+xr*xextend)
  yr <- c(min(y)-yr*yextend,max(y)+yr*yextend)
  if(kde.package=="ash")
  {
    if(is.null(h))
      h <- c(5,5)
    den <- ash::ash2(ash::bin2(cbind(x,y),rbind(xr,yr)),h)
  }
  else
  {
    X <- cbind(x,y)
    if(is.null(h))
      h <- ks::Hpi.diag(X,binned=TRUE)
    else
      h <- diag(h)
    den <- ks::kde(x=X,H=h,xmin=c(xr[1],yr[1]),xmax=c(xr[2],yr[2]))
    den <- list(x=den$eval.points[[1]],y=den$eval.points[[2]],z=den$estimate)
  }
  return(den)
}

