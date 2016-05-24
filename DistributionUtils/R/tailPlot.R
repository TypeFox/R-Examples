tailPlot <- function(x, log = "y", side = c("right", "left"),
                     main = NULL, xlab = NULL, ylab = NULL, ...)
{
  ## Purpose: Draws a tail plot
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 10 Mar 2010, 16:01

  xName <- paste(deparse(substitute(x), 500), collapse = "\n")
  side <- match.arg(side)
  x <- sort(as.numeric(x))
  ypoints <- 1 - ppoints(x)
  med <- median(x)
  xRange <- range(x)
  if (side =="left"){
    ypoints <- ppoints(x)
    xlim <- c(xRange[1],med)
    if (is.null(xlab)){
      xlab <- "x"
      if (log == "x" || log == "xy" || log == "xy"){
        xlab <- paste(xlab, "(on log scale)")
      }
    }
    if (is.null(ylab)){
      ylab <- "F(x)"
      if (log == "y" || log == "xy" || log == "xy"){
        ylab <- paste(ylab, "(on log scale)")
      }
    }
    if (is.null(main)){
      main = paste("Left Tail Plot of", xName)
    }
    plot(x, ypoints, log = log, xlim = xlim, main = main,
         xlab = xlab, ylab = ylab, ...)
  }
  if (side =="right"){
    ypoints <- 1 - ppoints(x)
    xlim <- c(med,xRange[2])
    if (is.null(xlab)){
      xlab <- "x"
      if (log == "x" || log == "xy" || log == "xy"){
        xlab <- paste(xlab, "(on log scale)")
      }
    }
    if (is.null(ylab)){
      ylab <- "F(x)"
      if (log == "y" || log == "xy" || log == "xy"){
        ylab <- paste(ylab, "(on log scale)")
      }
    }
    if (is.null(main)){
      main = paste("Right Tail Plot of", xName)
    }
    plot(x, ypoints, log = log, xlim = xlim, main = main,
         xlab = xlab, ylab = ylab, ...)
  }
  invisible(NULL)
}

#----------------------------------------------------------

tailPlotLine <- function(x, distrFn, param = NULL,
                         side = c("right", "left"), ...)
{
  ## Purpose: Add distribution line to tail plot for distributions
  ## which use the standardised approach
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 11 Mar 2010, 11:35
  
  pfun <- match.fun(paste("p", distrFn, sep = ""))
  side <- match.arg(side)
  x <- sort(as.numeric(x))
  if (side =="left"){
    if (!is.null(param)){
      lines(x, pfun(x, param = param), ...)
    } else {
      lines(x, pfun(x), ...)
    }
  }
  if (side =="right"){
    if (!is.null(param)){
      lines(x, 1 - pfun(x, param = param), ...)
    } else {
      lines(x, 1 - pfun(x), ...)
    }
  }
  invisible(NULL)
}


#----------------------------------------------------------

normTailPlotLine <- function(x, mean = 0, sd = 1,
                             side = c("right", "left"), ...)
{
  ## Purpose: Add normal distribution line to tail plot
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 10 Mar 2010, 15:58

  side <- match.arg(side)
  x <- sort(as.numeric(x))
  if (side =="left"){
    lines(x, pnorm(x, mean = mean, sd = sd), ...)
  }
  if (side =="right"){
    lines(x, 1 - pnorm(x, mean = mean, sd = sd), ...)
  }
  invisible(NULL)
}

#----------------------------------------------------------

tTailPlotLine <- function(x, df = Inf,
                          side = c("right", "left"), ...)
{
  ## Purpose: Add t distribution line to tail plot
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 10 Mar 2010, 15:58

  side <- match.arg(side)
  x <- sort(as.numeric(x))
  if (side =="left"){
    lines(x, pt(x, df = df), ...)
  }
  if (side =="right"){
    lines(x, 1 - pt(x, df = df), ...)
  }
  invisible(NULL)
}

#----------------------------------------------------------

gammaTailPlotLine <- function(x, shape = 1, rate = 1, scale = 1/rate,
                              side = c("right", "left"), ...)
{
  ## Purpose: Add gamma distribution line to tail plot
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## ----------------------------------------------------------------------
  ## Author: David Scott, Date: 10 Mar 2010, 15:58

  side <- match.arg(side)
  x <- sort(as.numeric(x))
  if (side =="left"){
    lines(x, pgamma(x, shape = shape, scale = scale), ...)
  }
  if (side =="right"){
    lines(x, 1 - pgamma(x, shape = shape, scale = scale), ...)
  }
  invisible(NULL)
}


