#' Plot conditional detection function from distance sampling model
#'
#' Plot proportion of observations detected within distance intervals (for
#' conditional detection functions) to compare visually the fitted model and
#' data. Internal function called by \code{plot} methods.
#'
#' @aliases plot_cond
#' @param obs obsever code
#' @param xmat processed data
#' @param gxvalues detection function values for each observation
#' @param model   fitted model from \code{ddf}
#' @param nc number of equal-width bins for histogram
#' @param breaks user define breakpoints
#' @param finebr fine break values over which line is averaged
#' @param showpoints logical variable; if \code{TRUE} plots predicted value
#'   for each observation
#' @param showlines logical variable; if \code{TRUE} plots average predicted
#'   value line
#' @param maintitle main title line for each plot
#' @param ylim range of y axis (default \code{c(0,1)})
#' @param angle shading angle for hatching
#' @param density shading density for hatching
#' @param col plotting colour
#' @param jitter scaling option for plotting points.  Jitter is applied to
#'   points by multiplying the fitted value by a random draw from a normal
#'   distribution with mean 1 and sd jitter.
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param subtitle if TRUE, shows plot type as sub-title
#' @param \dots other graphical parameters, passed to the plotting functions
#'   (\code{plot}, \code{hist}, \code{lines}, \code{points}, etc)
#' @return NULL
#' @author Jeff Laake, Jon Bishop, David Borchers
#' @keywords plot
plot_cond <- function(obs,xmat,gxvalues,model,nc,breaks,finebr,showpoints,
                      showlines,maintitle,ylim,angle=-45,density=20,col="black",
                      jitter=NULL,xlab="Distance",ylab="Detection probability",
                      subtitle=TRUE,...){

  # build and plot the histogram
  selection <- xmat$detected[xmat$observer!=obs]==1
  selmat <- (xmat[xmat$observer==obs,])[selection,]
  shist <- hist(xmat$distance[xmat$observer!= obs & xmat$detected==1],
                breaks=breaks, plot = FALSE)
  mhist <- hist(xmat$distance[xmat$timesdetected== 2 & xmat$observer==obs],
                breaks=breaks, plot = FALSE)
  if(length(mhist$counts) < length(shist$counts)){
    prop <- c(mhist$counts/shist$counts[1:length(mhist$counts)],
              rep(0, (length(shist$counts) - length(mhist$counts))))
  }else{
    prop <- mhist$counts/shist$counts
  }

  mhist$density <- prop
  mhist$equidist <- FALSE
  mhist$intensities <- mhist$density
  histline(mhist$density,breaks=breaks,lineonly=FALSE,xlab=xlab,ylab=ylab,
           ylim=ylim,fill=TRUE, angle=angle,density=density,col=col,
           det.plot=TRUE,...)

  # plot the detection function
  if(showlines){
    line <- average.line.cond(finebr,obs,model)
    linevalues <- line$values
    xgrid <- line$xgrid
    lines(xgrid,linevalues,...)
  }

  # plot points
  if(showpoints){
    ifelse(is.null(jitter),
           jitter.p<-1,
           jitter.p<-rnorm(length(gxvalues),1,jitter))
    points(selmat$distance,gxvalues*jitter.p,...)
  }

  # add the usual titles as subtitles if that's what the user asked for
  if(maintitle!=""){
    if(subtitle){
      title(paste(maintitle,
                  "\nConditional detection probability\nObserver=",
                  obs ," | Observer = ",3-obs),...)
    }else{
     title(maintitle)
    }
  }else{
    if(subtitle){
      title(paste("Conditional detection probability\nObserver=",
                  obs ," | Observer = ",3-obs),...)
    }
  }
}
