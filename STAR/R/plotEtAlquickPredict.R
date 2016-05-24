##############################################################################
##############################################################################
##############################################################################
plot.quickPredict <- function(x,
                              y,
                              xFct=function(x) x,
                              style=c("band",
                                "simple"),
                              ylim,
                              meanCol=2,
                              meanWD=2,
                              bandCol="grey50",
                              xlab,
                              ylab,
                              ...
                              )
{
#######################################################################
### Method / function plot.quickPredict
### ----------------------------------------------------------
### Arguments:
###  x: a "quickPredict" object.
###  y: not used, only included for compatibility with generic method.
###  xFct: a function applied to the x axis.
###  style: a character, either "band" or "simple". Controls the way
###         confidence intervals are displayed. If "band" is selected
###         they are shown as a ribbon. Their boundaries appear as
###         dashed lines otherwise.
### ylim: see plot (default provided if missing).
### meanCol: the color used to display the estimated mean of the term
###          (see argument "col" in plot).
### meanWD: the width used to display the estimated mean of the term
###         (see argument "lwd" in plot).
### bandCol: the color of the confidence interval ribbon when one is
###          is drawn.
### xlab, ylab: see plot (default provided if any is missing).
### ...:  see plot.
#######################################################################
  xx <- x$xx
  est.mean <- x$est.mean
  
  ylimMissing <- missing(ylim)
  if (ylimMissing) ylim <- range(est.mean)
  se.fit <- !is.null(x$est.sd)
  if (se.fit) {
    est.upr <- est.mean + 1.96 * x$est.sd
    est.lwr <- est.mean - 1.96 * x$est.sd
    if (ylimMissing) ylim <- c(min(est.lwr),max(est.upr))
  } ## End of conditional on se.fit

  xx <- xFct(xx)
  if (missing(xlab)) xlab <- x$include
  if (missing(ylab)) ylab <- expression(eta)
  plot(xx,
       est.mean,
       type="n",
       ylim=ylim,
       xlab=xlab,
       ylab=ylab,
       ...)
  if (se.fit) {
    if (style[1] == "simple") {
      lines(xx,est.mean,col=meanCol,lwd=meanWD)
      lines(xx,est.upr,lty=2)
      lines(xx,est.lwr,lty=2)
    } ## End of conditional on style[1] == "simple"
    if (style[1] != "simple") {
      polygon(c(xx,rev(xx)),
              c(est.upr,rev(est.lwr)),
              border=NA,
              col=bandCol)
      lines(xx,est.mean,col=meanCol,lwd=meanWD)
    } ## End of conditional on style[1] != "simple"
  } else {
    lines(xx,est.mean,col=meanCol,lwd=meanWD)
  } ## End of conditional on se.fit
  
}
##############################################################################
##############################################################################
##############################################################################
lines.quickPredict <- function(x,
                               what=c("mean","upr","lwr"),
                               xFct=function(x) x,
                               ...
                               )
{
#######################################################################
### Method / function lines.quickPredict
### ----------------------------------------------------------
### Arguments:
###  x: a "quickPredict" object.
###  what: one of the following character strings: "mean", "upr", "lwr".
###        Controls the line drawn, the estimated mean, upper bound of
###        the 95% CI or lower bound.
###  xFct: a function applied to the x axis.
### ...:  see lines.
#######################################################################
  if (what[1] != "mean" && is.null(x$est.sd))
    stop("Standard deviation of the estimate missing.\n Re-run quickPredict with \"se.fit\" set to TRUE.")
  
  xx <- x$xx
  
  xx <- xFct(xx)
  switch(what[1],
         mean = lines(xx,x$est.mean,...),
         upr = lines(xx,x$est.mean + 1.96 * x$est.sd,...),
         lwr = lines(xx,x$est.mean - 1.96 * x$est.sd,...)
         )
  
}
##############################################################################
##############################################################################
##############################################################################
image.quickPredict <- function(x,
                               main,
                               xlab,
                               ylab,
                               ...) {
#######################################################################
### Method / function image.quickPredict
### ----------------------------------------------------------
### Arguments:
###  x: a "quickPredict" object.
###  main, xlab, ylab: see image (default provided if any is missing).
###  ...: additional arguments of image.
#######################################################################
  ## Check that the quickPredict object x corresponds
  ## to an interaction term and that its yy component
  ## is not null
  if (is.null(x$yy)) stop("Component yy of x should be non NULL.")

  if (missing(main)) main <- paste("term",x$include)
  vNames <- strsplit(x$include,":")[[1]]
  if (missing(xlab)) xlab <- vNames[1]
  if (missing(ylab)) ylab <- vNames[2]

  zz <- x$est.mean
  xx <- x$xx
  yy <- x$yy
  
  image(xx,
        yy,
        zz,
        xlab=xlab,
        ylab=ylab,
        main=main,...)
  
}
##############################################################################
##############################################################################
##############################################################################
contour.quickPredict <- function(x,
                                 what=c("mean","sd"),
                                 main,
                                 xlab,
                                 ylab,
                                 add,
                                 ...) {
#######################################################################
### Method / function contour.quickPredict
### ----------------------------------------------------------
### Arguments:
###  x: a "quickPredict" object.
###  what: a character string specifying if the mean or the sd contours
###        should be drawn.
###  main, xlab, ylab, add: see contour (default provided if any is missing).
###  ...: additional arguments of contour.
#######################################################################
  ## Check that the quickPredict object x corresponds
  ## to an interaction term and that its yy component
  ## is not null
  if (is.null(x$yy)) stop("Component yy of x shoulc be non NULL.")

  vNames <- strsplit(x$include,":")[[1]]
  if (missing(add)) add <- FALSE
  if (missing(xlab) && !add) xlab <- vNames[1]
  if (missing(ylab) && !add) ylab <- vNames[2]
  if (missing(main) && !add) main <- paste("term",x$include)

  zz <- switch(what[1],
               mean=x$est.mean,
               sd=x$est.sd
               )
  xx <- x$xx
  yy <- x$yy
  
  if (!add) {
    contour(xx,
            yy,
            zz,
            xlab=xlab,
            ylab=ylab,
            main=main,
            add=add,
            ...)
  } else {
    contour(xx,
            yy,
            zz,
            add=add,
            ...)
  }
  
}
##############################################################################
##############################################################################
##############################################################################
persp.quickPredict <- function(x,
                               what=c("mean","sd"),
                               main,
                               xlab,
                               ylab,
                               zlab,
                               ...) {
#######################################################################
### Method / function persp.quickPredict
### ----------------------------------------------------------
### Arguments:
###  x: a "quickPredict" object.
###  what: a character string specifying if the mean or the sd perspective
###        should be drawn.
###  main, xlab, ylab, zlab: see persp (default provided if any is missing).
###  ...:  additional arguments of persp.
#######################################################################
  ## Check that the quickPredict object x corresponds
  ## to an interaction term and that its yy component
  ## is not null
  if (is.null(x$yy)) stop("Component yy of x shoulc be non NULL.")

  vNames <- strsplit(x$include,":")[[1]]
  if (missing(xlab)) xlab <- vNames[1]
  if (missing(ylab)) ylab <- vNames[2]
  if (what[1] == "mean") {
    if (missing(zlab)) zlab <- "mean"
    if (missing(main)) main <- paste("Mean of term",x$include)
  } else {
    if (missing(zlab)) zlab <- "sd"
    if (missing(main)) main <- paste("SD of term",x$include)
  }

  zz <- switch(what[1],
               mean=x$est.mean,
               sd=x$est.sd
               )
  xx <- x$xx
  yy <- x$yy
  
  persp(xx,
        yy,
        zz,
        xlab=xlab,
        ylab=ylab,
        zlab=zlab,
        main=main,
        ...)
  
}
