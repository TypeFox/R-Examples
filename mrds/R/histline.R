#' Plot histogram line
#'
#' Takes bar heights (height) and cutpoints (breaks), and constructs a
#' line-only histogram from them using the function plot() (if lineonly==FALSE)
#' or lines() (if lineonly==TRUE).
#'
#' @param height heights of histogram bars
#' @param breaks cutpoints for x
#' @param lineonly if TRUE, drawn with plot; otherwise with lines to allow
#'   addition of current plot
#' @param outline if TRUE, only outline of histogram is plotted
#' @param fill If fill==TRUE, uses polygon() to fill bars
#' @param ylim limits for y axis
#' @param xlab label for x axis
#' @param ylab label for y axis
#' @param det.plot if TRUE, plot is of detection so yaxis limited to unit
#'   interval
#' @param add should this plot add to a previous window
#' @param \dots Additional unspecified arguments for plot (fill==TRUE
#' @return None
#' @author ???
#' @importFrom graphics plot box
histline <- function(height, breaks, lineonly=FALSE, outline=FALSE, fill=FALSE,
                     ylim=range(height), xlab="x", ylab="y", det.plot=FALSE,
                     add=FALSE,...){

  # make hist object
  hh <- list()
  hh$breaks <- breaks
  hh$counts <- height
  hh$density <- height#/sum(height)
  hh$mids <- breaks[-length(breaks)] + diff(breaks)/2
  hh$xname <- "hh"
  hh$equidist <- FALSE
  class(hh) <- "histogram"

  # make the plot
  if(det.plot){
    plot(hh,ylim=ylim,xlab=xlab,ylab=ylab,yaxp=c(0,1,5),main="",add=add,freq=FALSE,...)
  }else{
    plot(hh,ylim=ylim,xlab=xlab,ylab=ylab,main="",add=add,freq=FALSE,...)
  }
  # put a box around
  box()
}
