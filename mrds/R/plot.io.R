#' Plot fit of detection functions and histograms of data from distance
#' sampling independent observer (\code{io}) model
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
#' @export
#' @param x fitted model from \code{ddf}
#' @param which index to specify which plots should be produced.
#'  \tabular{ll}{1 \tab Plot primary unconditional detection function \cr
#'               2 \tab Plot secondary unconditional detection function \cr
#'               3 \tab Plot pooled unconditional detection function \cr
#'               4 \tab Plot duplicate unconditional detection function \cr
#'               5 \tab Plot primary conditional detection function\cr
#'               6 \tab Plot secondary conditional detection function \cr}
#'  Note that the order of which is ignored and plots are produced in the above
#'  order.
#' @param breaks user define breakpoints
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
#' @return Just plots
#' @author Jeff Laake, Jon Bishop, David Borchers, David L Miller
#' @keywords plot
#' @examples
#' \donttest{
#' library(mrds)
#' data(book.tee.data)
#' egdata <- book.tee.data$book.tee.dataframe
#' result.io <- ddf(dsmodel=~cds(key = "hn"), mrmodel=~glm(~distance),
#'                  data=egdata, method="io", meta.data=list(width=4))
#'
#' # just plot everything
#' plot(result.io)
#'
#' # Plot primary and secondary unconditional detection functions on one page
#' # and  primary and secondary conditional detection functions on another
#' plot(result.io,which=c(1,2,5,6),pages=2)
#' }
plot.io <- function(x, which=1:6, breaks=NULL, nc=NULL,  maintitle="",
                    showlines=TRUE, showpoints=TRUE,ylim=c(0,1),angle=-45,
                    density=20,col="black",jitter=NULL,divisions=25,pages=0,
                    xlab="Distance",ylab="Detection probability",subtitle=TRUE,
                    ...){

  model <- x

  # since we can't get the same ordering as which is, make sure it's the same
  # every time, at least
  which <- sort(which)

  # Retrieve values from model object
  xmat.p0 <- model$mr$mr$data
  xmat.p0$offsetvalue <- 0
  xmat.p0$distance <- 0
  ddfobj <- model$ds$ds$aux$ddfobj
  if(ddfobj$type=="gamma"){
    xmat.p0$distance <- rep(apex.gamma(ddfobj),2)
  }

  p0 <- predict(model$mr,newdata=xmat.p0,integrate=FALSE)$fitted
  xmat <- model$mr$mr$data
  cond.det <- predict(model$mr,newdata=xmat,integrate=FALSE)
  width <- model$meta.data$width
  left <- model$meta.data$left
  detfct.pooled.values <- detfct(xmat$distance[xmat$observer==1],
                                 ddfobj,width=width-left)
  delta <- cond.det$fitted/(p0*detfct.pooled.values)
  p1 <- cond.det$p1
  p2 <- cond.det$p2

  # list of values to pass as gxvalues for unconditional plots
  # in order as below
  # gxvalues are the points shown in the plots
  gxlist <- list(p1/delta,
                 p2/delta,
                 (p1+p2-p1*p2)/delta,
                 p1*p2/delta)

  # If number of classes for histogram intervals was not set compute
  # a reasonable default
  if(is.null(nc)){
    nc<-round(sqrt(min(length(xmat$distance[xmat$observer==1&xmat$detected==1]),
                       length(xmat$distance[xmat$observer==2&xmat$detected==1]),
                       length(xmat$distance[xmat$observer==1 &
                              xmat$timesdetected==2]) )),0)
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

  # loop over the unconditional plots
  # 1 - Plot primary unconditional detection function
  # 2 - Plot secondary unconditional detection function
  # 3 - Plot pooled unconditional detection function
  # 4 - Plot duplicate unconditional detection function
  for(wh in which[which<5]){
    plot_uncond(model, wh, xmat, gxvalues=gxlist[[wh]], nc,
                finebr=(width/divisions)*(0:divisions), breaks, showpoints,
                showlines, maintitle, ylim,
                angle=angle,density=density,col=col,jitter=jitter,xlab=xlab,
                ylab=ylab,subtitle=subtitle,...)
  }

  # 5 - Plot conditional detection function (1|2)
  data <- model$mr$mr$data
  data$offsetvalue <- 0
  if(is.element(5,which)){
    gxvalues <-p1[xmat$detected[xmat$observer==2]==1]
    plot_cond(1,data,gxvalues,model,nc,breaks,
              finebr=(width/divisions)*(0:divisions),showpoints,showlines,
              maintitle,ylim,angle=angle,density=density,col=col,jitter=jitter,
              xlab=xlab,ylab=ylab,subtitle=subtitle,...)
  }

  # 6 - Plot secondary conditional detection function (2|1)
  if(is.element(6,which)){
    gxvalues <- p2[xmat$detected[xmat$observer==1]==1]
    plot_cond(2,data,gxvalues,model,nc,breaks,
              finebr=(width/divisions)*(0:divisions),showpoints,showlines,
              maintitle,ylim,angle=angle,density=density,col=col,jitter=jitter,
              xlab=xlab,ylab=ylab,subtitle=subtitle,...)
  }
  invisible(NULL)
}
