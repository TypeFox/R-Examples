#-------------------------------------------------------------------------
# plot.mqcs
#-------------------------------------------------------------------------
##' Plot method for 'mqcs' objects
##' 
##' Generic function for plotting Multivarite charts of object of class 'mqcs' to perform statistical 
##' quality control.
##' 
##' @method plot mqcs
##' @param x  Object mqcs (Multivarite Quality Control Statical)
##' @param title an overall title for the plot
##' @param subtitle a sub title for the plot
##' @param xlab a title for the x axis
##' @param ylab a title for the y axis
##' @param ylim the y limits of the plot
##' @param ...  arguments to be passed to or from methods.
##' @export
##' 
plot.mqcs <- function(x, title, subtitle, xlab, ylab, ylim, ...)
  #.........................................................................                     
{
  
  oldpar <- par(mar = c(5, 4, 4, 3) + 0.1)
  limits <- x$limits
  sample <- 1:length(x$statistics)

  plot(x$statistics ~ sample, type =  "n",pch  =  16, axes  =  FALSE, 
         main  =  title, sub  =  subtitle, xlab  = xlab, 
         ylab  =  ylab, ylim  =  ylim)    
  
  axis(1, at  =  sample, cex.axis  =  0.7)
  axis(2, cex.axis  =  0.7)
  
  axis(4, at = c(max(x$limits), min(x$limits)), 
       labels = c("UCL","LCL"), adj = 0, las = 1)
    
  
  rect(par("usr")[1],
       par("usr")[3],
       par("usr")[2],
       par("usr")[4],
       col  =  "#CCCCCC")
  box(col  =  "#CCCCCC")
  grid(col  =  "#EEEEEE")
  
  
      points(sample, x$statistics, pch = 3, cex = 0.8)
      lines(x$statistics ~ sample, type = "o", pch=20)      


    lcl <- x$limits[1]
    ucl <- x$limits[2]
    
    abline(h  =  lcl, lwd  =  2, col  =  "red", lty = 2)
    abline(h  =  ucl, lwd  =  2, col  =  "red",lty = 2) 
      
  
    beyond.limits <- x$violations

    
    points(x$statistics[beyond.limits]~sample[beyond.limits],
           col = "red",
           pch  =  19)
    
  par(oldpar)
  #.........................................................................
} # plot.mqcs
#.........................................................................

#-------------------------------------------------------------------------
# plot.mqcs.t2
#-------------------------------------------------------------------------
##' @rdname  plot.mqcs
##' @method  plot mqcs.t2
##' @export
##' 

plot.mqcs.t2 <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                         ylab  =  NULL, ylim  =  NULL, ...)
  #.........................................................................                     
{
  
    
  if(is.null(ylim)) 
    ylim <-  range(x$statistics, x$limits)
  
  if (is.null(title)) title <- expression(paste("Chart of control ", T ^ 2," "))
  
  if (is.null(subtitle)) subtitle <- "Hotelling"
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- expression(T ^ 2)
  
  plot.mqcs(x, title , subtitle, xlab, ylab, ylim)
} #plot.mqcs.t2
#.........................................................................                     


#-------------------------------------------------------------------------
# plot.mqcs.mcusum
#-------------------------------------------------------------------------
##' @rdname  plot.mqcs
##' @method  plot mqcs.mcusum
##' @export
##' 
plot.mqcs.mcusum <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                            ylab  =  NULL, ylim  =  NULL, ...)
  #.........................................................................                     
{
  
  
  if(is.null(ylim)) 
    ylim <-  range(x$statistics, x$limits)
  
  if (is.null(title)) title <- "Chart of control MCUSUM"
  
  if (is.null(subtitle)) subtitle <- "MCUSUM"
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- expression(mcusum)
  
  plot.mqcs(x, title , subtitle, xlab, ylab, ylim)
} #plot.mqcs.mcusum
#.........................................................................                     

#-------------------------------------------------------------------------
# plot.mqcs.mewma
#-------------------------------------------------------------------------
##' @rdname  plot.mqcs
##' @method  plot mqcs.mewma
##' @export
##' 

plot.mqcs.mewma <- function(x, title  =  NULL, subtitle  =  NULL, xlab  =  NULL,
                         ylab  =  NULL, ylim  =  NULL, ...)
  #.........................................................................                     
{
  
  
  if(is.null(ylim)) 
    ylim <-  range(x$statistics, x$limits)
  
  if (is.null(title)) title <- "Chart of control MEWMA"
  
  if (is.null(subtitle)) subtitle <- "MEWMA"
  if (is.null(xlab)) xlab <- "Sample"
  if (is.null(ylab)) ylab <- expression(mewma)
  
  plot.mqcs(x, title , subtitle, xlab, ylab, ylim)
} #plot.mqcs.mewma
#.........................................................................                     