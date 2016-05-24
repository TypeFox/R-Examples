scatterplot.density <- function(x, y, zlim, xylim, num.bins=64,
                                col=kristen.colors(32), xlab, ylab, main,
                                density.in.percent=TRUE,
                                col.regression.line=1,
                                col.one.to.one.line=grey(0.4),
                                col.bar.legend=TRUE, 
                                plt.beyond.zlim=FALSE, ...){

  ## Both and x and y should be numeric.
  if (!is.numeric(x))
    stop("x must be numeric")
  if (!is.numeric(y))
    stop("y must be numeric")



  ## x and y can be in a matrix or data.frame, but these will be forced
  ## into vector format.  The lengths of these vectors must be the same.
  x.data <- as.vector(x)
  y.data <- as.vector(y)
  if (length(x.data)!=length(y.data))
    stop("x and y must have the same length.")



  ## Determine the number of bins.
  if (is.list(num.bins)){  
    num.bins.x <- num.bins$num.bins.x
    num.bins.y <- num.bins$num.bins.y   
  }
  else{
    num.bins.x <- num.bins
    num.bins.y <- num.bins
  }



  ## Determine the x and y limits of the plot.
  if (missing(xylim)){ 
    ## When limits are missing, select limits which include all the data.
    xlim <- range(x.data)
    ylim <- range(y.data)
  }
  else{
    ## Use given limits.
    if (is.list(xylim)) {
      xlim <- xylim$xlim
      ylim <- xylim$ylim
    }
    else{
      xlim <- xylim
      ylim <- xylim
    }
  }



  ## Find suitable bins, making sure to enclose the xlim and ylim
  ## values in bins.
  data.bins.x <- seq(xlim[1], xlim[2], length=num.bins.x)
  bin.x.length <- data.bins.x[2] - data.bins.x[1]
  plot.seq.x <- seq(xlim[1]-(bin.x.length/2), xlim[2]+(bin.x.length/2), length=num.bins.x+1)

  data.bins.y <- seq(ylim[1], ylim[2], length=num.bins.y)
  bin.y.length <- data.bins.y[2] - data.bins.y[1]
  plot.seq.y <- seq(ylim[1]-(bin.y.length/2), ylim[2]+(bin.y.length/2), length=num.bins.y+1)


  ## Figure out which bin each data point falls in.
  x.cut <- cut(x.data, plot.seq.x)
  y.cut <- cut(y.data, plot.seq.y)
  tab.x.y <- table(x.cut, y.cut)
  if (density.in.percent) 
    tab.x.y <- tab.x.y/length(x.data)*100
  ## If no data points fall in the bin, set that bin equal to NA so
  ## it doesn't get plotted with a color at all.
  tab.x.y[tab.x.y==0] <- NA



  ## If the user doesn't specify xlab and ylab, assign to them the
  ## expressions the user named as parameters x and y.  This is
  ## similar to the default behavior for image().
  if (missing(xlab))
    xlab <- deparse(substitute(x))
  if (missing(ylab))
    ylab <- deparse(substitute(y))

  ## If the user doesn't specify a main title, add one (based on
  ## whether the user has chosen to view densities in percent or in
  ## counts).
  if (missing(main)){
    if (density.in.percent)
      main <- "Data Density Plot (%)"
    else 
      main <- "Data Frequency Plot (counts)"
  }



  ## If plt.beyond.zlim=TRUE and zlim is given, we set z values
  ## which are greater than zlim[2] equal to zlim[2] and z values
  ## which are less than zlim[1] equal to zlim[1].  This ensures that
  ## these values are plotted, when they would otherwise be blank.  If
  ## plt.beyond.zlim=TRUE and zlim is missing, then zlim will be
  ## determined according to the total range of the z's and
  ## plt.beyond.zlim can be set to FALSE.
  if (plt.beyond.zlim){
    if (missing(zlim)){
      warning("plt.beyond.zlim=TRUE is not a valid option if zlim argument is not provided, changing to plt.beyond.zlim=FALSE")
      plt.beyond.zlim <- FALSE
    }
    else{
      ## Values lower than zlim[1] set to zlim[1].
      tab.x.y[tab.x.y < zlim[1]] <- zlim[1]
      ## Values higher than zlim[2] set to zlim[2].
      tab.x.y[tab.x.y > zlim[2]] <- zlim[2]
    }
  }
  ## By default, zlim is set to minimum and maximum values in
  ## tab.x.y.  That is, we set zlim to the minimum and maximum
  ## number/percentage of data points in the cells.
  if (missing(zlim)) 
    zlim <- range(tab.x.y, na.rm=T)



  ## Draw color-coded square bins using image().
  image(x=plot.seq.x, y=plot.seq.y, z=tab.x.y, zlim=zlim, col=col,
        xlab=xlab, ylab=ylab, main=main, ...)



  ## If col.one.to.one.line flag is not null, include this line (slope 1 and
  ## intercept 0) with the color specified.
  if (!is.null(col.one.to.one.line))
    abline(0, 1, col=col.one.to.one.line, lty=3)

  ## If col.regression.line flag is not null, include the regression
  ## line with the color specified.
  if (!is.null(col.regression.line)){
    y.x.lm <- lm(y.data~x.data)$coeff
    abline(y.x.lm, col=col.regression.line)
    ## Build legend text that contains the regression equation.
    legend.txt <- bquote( paste( hat(y), "=", .(signif(y.x.lm[1], 2)),
                                "+", .(signif(y.x.lm[2], 2)), "x") )
    legend("topleft", legend=do.call("expression", c(legend.txt)), bty="n",
           text.col=col.regression.line)
  }



  ## Add vertically-oriented color bar to explain the correspondence
  ## between colors and percentages/counts.
  if (col.bar.legend)
    vertical.image.legend(col=col, zlim=zlim)
}
