conf2d <- function(x, ...)
{
  UseMethod("conf2d")
}


conf2d.formula <- function(formula, data, subset, ...)
{
  m <- match.call(expand.dots=FALSE)
  if(is.matrix(eval(m$data,parent.frame())))
    m$data <- as.data.frame(data)
  m$... <- NULL
  m[[1L]] <- as.name("model.frame")
  mf <- eval(m, parent.frame())

  conf2d.default(mf[2:1], ...)
}


conf2d.default <- function(x, y, level=0.95, n=200, method="wand", shape=1, smooth=50, plot=TRUE, add=FALSE, xlab=NULL,
                           ylab=NULL, col.points="gray", col="black", lwd=2, ...)
{
  method <- match.arg(tolower(method), c("wand","mass"))

  ## 1  Extract data
  if(is.matrix(x))
    x <- as.data.frame(x)
  if(is.list(x))  # data.frame or list
  {
    xlab <- if(is.null(xlab)) names(x)[1] else xlab
    ylab <- if(is.null(ylab)) names(x)[2] else ylab
    y <- x[[2]]
    x <- x[[1]]
  }

  ## 2  Plot xy scatter
  if(plot && !add)
  {
    if(is.null(xlab))
      xlab <- deparse(substitute(x))
    if(is.null(ylab))
      ylab <- deparse(substitute(y))
    plot(x, y, xlab=xlab, ylab=ylab, lwd=1, col=col.points, ...)
  }

  ## 3  Compute surface matrix
  smooth <- rep(smooth, length.out=2)
  if(method == "wand")
  {
    bandwidth <- c(dpik(x), dpik(y))
    surf <- bkde2D(cbind(x,y), bandwidth=1.3*bandwidth*shape, gridsize=c(smooth,smooth))
  }
  else  # mass
  {
    bandwidth <- c(bandwidth.nrd(x), bandwidth.nrd(y))
    lx <- range(x) + c(-0.05,0.05)*diff(range(x))
    ly <- range(y) + c(-0.05,0.05)*diff(range(y))
    surf <- kde2d(x, y, h=bandwidth*shape, n=5*smooth/shape, lims=c(lx,ly))
  }

  ## 4  Find best region
  output <- conf2d_int(x, y, surf, level, n)
  if(abs(level-output$prop) > 0.01)
     warning("appropriate region not found (desired level=", level, ", but actual prop=", round(output$prop,3), ")")

  ## 5  Overlay polygon
  if(plot || add)
  {
    polygon(output, lwd=lwd, border=col, col=NA, ...)
    return(invisible(output))
  }
  else
  {
    return(output)
  }
}


conf2d_int <- function(x, y, surf, level, n)
{
  ## 1  Slice surface into candidate regions
  cl <- contourLines(surf[[1]], surf[[2]], surf[[3]], nlevels=n)  # candidate regions as contour lines

  ## 2  Find best region by counting the points inside
  pts <- SpatialPoints(cbind(x,y))
  spols <- list()             # spatial polygons representing candidate regions
  pin <- numeric(length(cl))  # number of points inside each region
  for(i in seq_len(length(cl)))
  {
    spol <- tryCatch(Polygon(cbind(cl[[i]]$x,cl[[i]]$y)), error=function(...) NA)
    spol <- tryCatch(Polygons(list(spol),ID=" "), error=function(...) NA)
    spol <- tryCatch(SpatialPolygons(list(spol)), error=function(...) NA)
    pin[i] <- tryCatch(sum(!is.na(over(pts,spol))), error=function(...) 0)
    spols[[i]] <- spol
  }
  pin <- pin / length(x)
  best <- which.min(abs(pin-level))

  ## 3  Extract statistics
  best.spol <- spols[[best]]@polygons[[1]]@Polygons[[1]]
  xcoords <- best.spol@coords[,1]
  ycoords <- best.spol@coords[,2]
  inside <- !is.na(over(pts,spols[[best]]))
  area <- best.spol@area
  prop <- pin[best]
  output <- list(x=xcoords, y=ycoords, inside=inside, area=area, prop=prop)

  return(output)
}
