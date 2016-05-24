#' Root Square Difference
#' @param x a column vector of scores on which the rsd is conditioned
#' @param o a column vector of equated scores based on the overall population (aligned with elements in x)
#' @param g a column vector of equated scores based on a single subpopulation (aligned with elements in x)
#' @param d a scalar of the difference that matters
#' @param s a scalar representing the standard deviation of x for any (sub)population of interest (e.g., synthetic population) (default is 1, which leads to calculation of the unstandardized rsd)
#' @param ymax A maximum value for the y axis of the plot (default = 3 + the maximum RSD value)
#' @param xlab A label for the x axis of the plot (default = Score Scale)
#' @param color of RSD line and points on plot (default = red)
#' @author Anne Corinne Huggins-Manley 
#' @references 
#' \itemize{
#' \item{Huggins, A.C., & Penfield, R.D. (2012). An NCME instructional module on population invariance in linking and equating. Educational Measurement: Issues and Practices, 31, 27-40.}
#' \item{Yang, W.L. (2004). Sensitivity of linkings between AP multiple-choice scores and composite scores to geographical region: An illustration of checking for population invariance. Journal of Educational Measurement, 41, 33-41.}
#' }
#' @seealso \code{\link{resd}}
#' @examples
#' #Unstandardized RSD for subpopulation 1 in the example data set, ex.data
#' rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,3],d=.5)
#' 
#' #Unstandardized RSD for subpopulation 5 in the example data set, ex.data
#' rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],d=.5)
#' 
#' #Unstandardized RSD for subpopulation 5 in the example data set, ex.data with adjustments to the maximum y-axis on the plot, a new x label, and points/line in green
#' rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],d=.5,ymax=3,xlab="Old Form",color="green")
#' 
#' #Standardized RSD for subpopulation 5 in the example data set, ex.data
#' rsd(x=ex.data[,1],o=ex.data[,2],g=ex.data[,7],d=.5,s=4.2)
#' @return A data frame of root square difference indices, conditioned on the score scale
#' @return A plot of the RSD in reference to the difference that matters
#' @export

rsd<- function(x,o,g,d,s,ymax,xlab,color){
  if(missing(s))
    s <- 1
  
  rsd<-(sqrt((o-g)^2))/s
  dfrsd <- data.frame("x"=x,"rsd"=rsd)
  
  if(missing(ymax))
    ymax <- max((dfrsd[,2])+3)
  if(missing(xlab))
    xlab <- "Score Scale"
  if(missing(color))
    color <- "red"
  
  plot<- plot(type="b",x=dfrsd[,1],y=dfrsd[,2],ylim=c(0,ymax),xlab=xlab,ylab="RSD",col=color,pch=19)
  abline(d/s,0,col="black",lty=2)
  legend(x=min(dfrsd[,1]),y=ymax,c("RSD","DTM"),col=c(color,"black"),pch=c(19,-1),lty=c(0,5))
  
  return(dfrsd)
  return(plot)
  
}