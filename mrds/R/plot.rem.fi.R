#' Plot fit of detection functions and histograms of data from removal distance
#' sampling model
#'
#' Plots the fitted detection functions for a distance sampling model and
#' histograms of the distances (for unconditional detection functions) or
#' proportion of observations detected within distance intervals (for
#' conditional detection functions) to compare visually the fitted model and
#' data.
#'
#' The structure of the histogram can be controlled by the user-defined
#' arguments \code{nc} or \code{breaks}.  The observation specific detection
#' probabilities along with the line representing the fitted average detection
#' probability.
#'
#' It is not intended for the user to call \code{plot.rem.fi} but its arguments
#' are documented here. Instead the generic \code{plot} command should be used
#' and it will call the appropriate function based on the class of the
#' \code{ddf} object.
#'
#' @aliases plot.rem.fi
#' @export
#' @param x fitted model from \code{ddf}
#' @param which index to specify which plots should be produced.
#'  \tabular{ll}{1 \tab Plot primary unconditional detection function \cr
#'               2 \tab Plot pooled unconditional detection function \cr
#'               3 \tab Plot conditional (1|2) detection function \cr}
#' @param breaks user defined breakpoints
#' @param nc number of equal-width bins for histogram
#' @param maintitle main title line for each plot
#' @param showpoints logical variable; if TRUE plots predicted value for each
#'   observation
#' @param showlines logical variable; if TRUE a line representing the average
#'   detection probability is plotted
#' @param ylim range of y axis; defaults to (0,1)
#' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col plotting colour
#' @param jitter scaling option for plotting points.  Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter
#' @param divisions number of divisions for averaging line values; default = 25
#' @param pages the number of pages over which to spread the plots. For
#'  example, if \code{pages=1} then all plots will be displayed on one page.
#'  Default is 0, which prompts the user for the next plot to be displayed.
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param subtitle if TRUE, shows plot type as sub-title
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (\code{plot}, \code{hist}, \code{lines}, \code{points}, etc)
#' @return NULL
#' @author Jeff Laake, Jon Bishop, David Borchers, David L Miller
#' @keywords plot
plot.rem.fi <- function(x,which=1:3,breaks=NULL,nc=NULL,maintitle="",
                        showlines=TRUE, showpoints=TRUE,ylim=c(0,1),angle=-45,
                        density=20,col="black",jitter=NULL,divisions=25,
                        pages=0,xlab="Distance",ylab="Detection probability",
                        subtitle=TRUE,...){
  # Functions used: process.data, predict(predict.io.fi),
  #                 plot.uncond, plot.cond

  # Retrieve values from model object
  model <- x
  xmat <- model$data
  xmat$offsetvalue <- 0
  cond.det <- predict(model,newdata=xmat,integrate=FALSE)
  fitted <- cond.det$fitted
  p1 <- cond.det$p1
  p2 <- cond.det$p2
  width <- model$meta.data$width
  left <- model$meta.data$left

  # If number of classes for histogram intervals was not set compute
  # a reasonable default
  if(is.null(nc)){
    nc <- round(sqrt(length(xmat$distance[xmat$observer==2&xmat$detected==1])))
  }

  # Set up default break points unless specified
  if(model$meta.data$binned){
    breaks <- model$meta.data$breaks
    nc <- length(breaks)-1
  }else{
    if(is.null(breaks)){
      breaks <- left + ((width-left)/nc)*(0:nc)
    }else{
      nc <- length(breaks)-1
    }
  }

  # do the plotting layout
  oask <- plot.layout(which,pages)
  on.exit(devAskNewPage(oask))

  # Plot primary unconditional detection function
  if(is.element(1,which)){
    plot_uncond(model,1,xmat,gxvalues=p1,nc,
                finebr=(width/divisions)*(0:divisions),breaks,showpoints,
                showlines,maintitle,ylim,
                angle=angle,density=density,col=col,jitter=jitter,
                xlab=xlab,ylab=ylab,subtitle=subtitle,...)
  }

  # Plot pooled unconditional detection function
  if(is.element(2,which)){
    plot_uncond(model,3,xmat,gxvalues=p1+p2*(1-p1),nc,
                finebr=(width/divisions)*(0:divisions),breaks,showpoints,
                showlines,maintitle,ylim,
                angle=angle,density=density,col=col,jitter=jitter,
                xlab=xlab,ylab=ylab,subtitle=subtitle,...)
  }

  # Plot conditional detection function
  data <- process.data(model$data,model$meta.data)$xmat
  data$offsetvalue <- 0
  if(is.element(3,which)){
    gxvalues <- p1[xmat$detected[xmat$observer==2]==1]
    plot_cond(1,data,gxvalues,model,nc,breaks,
              finebr=(width/divisions)*(0:divisions),showpoints,showlines,
              maintitle,ylim,angle=angle,density=density,col=col,jitter=jitter,
              xlab=xlab,ylab=ylab,subtitle=subtitle,...)
  }
  invisible(NULL)
}
