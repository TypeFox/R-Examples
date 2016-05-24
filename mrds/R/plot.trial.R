#' Plot fit of detection functions and histograms of data from distance
#' sampling trial observer model
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
#' It is not intended for the user to call \code{plot.io.fi} but its arguments
#' are documented here. Instead the generic \code{plot} command should be used
#' and it will call the appropriate function based on the class of the
#' \code{ddf} object.
#'
#' @aliases plot.trial
#' @export
#' @param x fitted model from \code{ddf}
#' @param which index to specify which plots should be produced.
#'  \tabular{ll}{1 \tab Unconditional detection function for observer 1 \cr
#'               2 \tab Conditional detection function plot (1|2)\cr}
#' @param breaks user define breakpoints
#' @param nc number of equal-width bins for histogram
#' @param maintitle main title line for each plot
#' @param showlines logical variable; if TRUE a line representing the average
#'   detection probability is plotted
#' @param showpoints logical variable; if TRUE plots predicted value for each
#'   observation
#' @param ylim range of y axis; defaults to (0,1)
#' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col plotting colour
#' @param jitter scaling option for plotting points.  Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter.
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
#' @author Jeff Laake, Jon Bishop, David Borchers
#' @keywords plot
plot.trial <- function(x, which=1:2, breaks=NULL, nc=NULL, maintitle="",
                       showlines=TRUE,showpoints=TRUE, ylim=c(0,1),angle=-45,
                       density=20,col="black",jitter=NULL,divisions=25,pages=0,
                       xlab="Distance",ylab="Detection probability",
                       subtitle=TRUE,...){
  # Uses: detfct, plot_cond, plot_uncond

  model <- x
  data <- process.data(model$data,model$meta.data)$xmat
  xmat <- data[data$observer==1&data$detected==1,]
  width <- model$meta.data$width
  left <- model$meta.data$left

  #  Set up default break points unless specified
  if(is.null(nc)){
    nc<-round(sqrt(min(length(xmat$distance[xmat$observer==1&xmat$detected==1]),
                       length(xmat$distance[xmat$observer==1 &
                              xmat$timesdetected==2]) )),0)
  }
  if(model$ds$meta.data$binned){
    breaks <- model$ds$ds$aux$breaks
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

  # Uncond detection fct for observer 1
  if(is.element(1,which)){
    ddfobj <- model$ds$ds$aux$ddfobj
    detfct.pooled.values <- detfct(xmat$distance,ddfobj,width=width-left)
    xmat.trial <- xmat
    xmat.trial$distance <- 0
    if(ddfobj$type=="gamma"){
      xmat.trial$distance <- as.vector(apex.gamma(ddfobj))
    }
    g0 <- predict(model$mr$mr,newdata=xmat.trial,type="response")
    gxvalues <- detfct.pooled.values*g0
    plot_uncond(model,1,xmat,gxvalues,nc,
                finebr=(width/divisions)*(0:divisions),breaks,showpoints,
                showlines,maintitle,ylim,
                angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,
                ylab=ylab,subtitle=subtitle,...)
  }

  # Conditional detection function plot: Obs 1 given obs 2.
  # Duplicate detections with MR model fitted and the estimated detection
  # probability from the MR part of the model
  if(is.element(2,which)){
    xmat <- data[data$observer==2&data$detected==1,]
    gxvalues <- predict(model$mr,newdata=xmat,type="response",
                        integrate=FALSE)$fitted

    plot_cond(1,data,gxvalues,model,nc,breaks,
              finebr=(width/divisions)*(0:divisions),showpoints,showlines,
              maintitle,ylim,angle=angle,density=density,col=col,jitter=jitter,
              xlab=xlab,ylab=ylab,subtitle=subtitle,...)

  }
  invisible(NULL)
}
