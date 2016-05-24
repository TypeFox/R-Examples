######################################################
## FRACTAL multidimensional kernel density estimator
##
## Class: KDE
## Constructor function: KDE
## Methods:
##
##   eda.plot.KDE
##   plot.KDE
##   print.KDE
##
######################################################

###
# KDE
###

"KDE" <- function(x, at=NULL, n.grid=100)
{
  data.name <- deparseText(substitute(x))

  # convert input to a matrix
  x <- as.matrix(x)

  # if at are NULL select
  # a multidimensional uniform grid
  # spanning the data
  if (is.null(at)){

    ndim   <- numCols(x)
    ranges <- apply(x, MARGIN=2, range)
    spans  <- vector(length=ndim, mode="list")

    for (i in seq(ndim)){
      spans[[i]] <- seq(from=ranges[1,i], to=ranges[2,i], length=n.grid)
    }

    at <- as.matrix(expand.grid(spans))
  }

  if (!is.matrix(at))
    stop("at must be a matrix")

  checkScalarType(n.grid,"integer")
  if (n.grid < 1)
    stop("n.grid must be positive")

  # estimate the KDE
  z <- as.vector(itCall("RS_fractal_kernel_density_estimate",
    x, at))
    #COPY=rep(FALSE,2),
    #CLASSES = rep("matrix", 2),
    #PACKAGE="ifultools"))

  oldClass(z) <- "KDE"
  attr(z, "training")  <- x
  attr(z, "at")        <- at
  attr(z, "data.name") <- data.name

  z
}

###
# eda.plot.KDE
###

"eda.plot.KDE" <- function(x, box=TRUE, adj.main=1, nlevels=10, labex=1, ...)
{
  # specify plot grid
  old.plt <- splitplot(2,2,1)
  on.exit(par(old.plt))

  if (numCols(attr(x,"training")) == 1){
    plot(x, ...)
    return(invisible(NULL))
  }
  xatt <- attributes(x)

  plot(x,style="original",add=TRUE)
  title("Training Data/Test Points", adj=adj.main)

  splitplot(2,2,2)
  plot(x,style="perspective",box=box,yaxp=c(0,max(x),1),add=TRUE)
  title("KDE Perspective", adj=adj.main)

  splitplot(2,2,3)
  plot(x,style="contour",labcex=par("cex")*labex,nlevels=nlevels,add=TRUE)
  title("KDE Contour", adj=adj.main)

  invisible(NULL)
}

###
# plot.KDE
###

"plot.KDE" <- function(x, style="original",
  dimensions=1:2, xlab=NULL, ylab=NULL, zlab="KDE", grid=FALSE, theta=120, phi=30, add=FALSE, ...)
{

  if (numCols(attr(x,"training")) == 1){
  	plot(as.vector(attr(x,"at")), asVector(x),
  	  xlab=attr(x,"data.name"), ylab="KDE", ...)
    return(invisible(NULL))
  }

  if (!add){
    old.plt <- splitplot(1,1,1)
    on.exit(par(old.plt))
  }
  old.xpd <- par(xpd=TRUE)
  on.exit(par(old.xpd))

  xatt     <- attributes(x)
  training <- xatt$training[,dimensions]
  at       <- xatt$at[,dimensions]
  ranges   <- apply(rbind(training, at), MARGIN=2, range)

  nms <- dimnames(training)[[2]]
  if (is.null(nms))
    nms <- c("X","Y")
  if (is.null(xlab))
    xlab <- nms[1]
  if (is.null(ylab))
    ylab <- nms[2]

  style <- match.arg(lowerCase(style),c("original","perspective","contour"))

  if (style == "original"){

    plot(training, type="p", xlim=ranges[,1], ylim=ranges[,2],
      xlab=xlab, ylab=ylab, ...)
    if (grid) points(at, pch=".")

  }
  else if (style == "perspective"){

      persp(akima::interp(at[,1], at[,2], asVector(x)),
        xlab=xlab, ylab=ylab, zlab=zlab,
        axes=TRUE, theta=theta, phi=phi,...)
  }
  else if (style == "contour"){

    contour(akima::interp(at[,1], at[,2], asVector(x)),
      xlab=xlab, ylab=ylab, ...)
  }

  invisible(NULL)
}

###
# print.KDE
###

"print.KDE" <- function(x, justify="left", sep=":", ...)
{
  xatt <- attributes(x)

  main <- paste("Kernel Density Function Estimate of", xatt$data.name, "data")

  z <- list(
    "Number of variables/dimensions"=numCols(xatt$at),
    "Number of evaluation points"=numRows(xatt$at),
    "Kernel"="Epanechnikov")

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

