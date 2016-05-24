#' Absolute Difference for Pairs - AD(x)
#' @param x a column vector of scores on which the AD(x) is conditioned
#' @param g1 a column vector of equated scores based on a single subpopulation (aligned with elements in x)
#' @param g2 a column vector of equated scores based on a different single subpopulation (aligned with elements in x)
#' @param d a scalar of the difference that matters 
#' @param s a scalar representing the standard deviation of x for any (sub)population of interest (e.g., synthetic population) (default is 1, which leads to calculation of the unstandardized adx)
#' @param ymax A maximum value for the y axis of the plot (default = 3 + the maximum AD(x) value)
#' @param xlab A label for the x axis of the plot (default = Score Scale)
#' @param color color of AD(x) line and points on plot (default = red)
#' @author Anne Corinne Huggins-Manley 
#' @references 
#' \itemize{
#' \item{Huggins, A.C., & Penfield, R.D. (2012). An NCME instructional module on population invariance in linking and equating. Educational Measurement: Issues and Practices, 31, 27-40.}
#' \item{Kolen, M.J., & Brennan, R.L. (2004). Test equating, scaling, and linking: Methods and practices (2nd ed.). NY: Springer.}
#' }
#' @seealso \code{\link{madp}}
#' @examples
#' #Unstandardized AD(x) for subpopulation 1 and subpopulation 2 in the example data set, ex.data
#' adx(x=ex.data[,1],g1=ex.data[,3],g2=ex.data[,4],d=.5)
#' 
#' #Unstandardized AD(x) for subpopulation 1 and subpopulation 2 in the example data set, ex.data, with adjustments to the maximum y-axis on the plot, a new xlabel, and points/line in blue.
#' adx(x=ex.data[,1],g1=ex.data[,3],g2=ex.data[,4],d=.5,ymax=2,xlab="AD(x)",color="blue")
#' 
#' #Unstandardized AD(x) for subpopulation 4 and subpopulation 5 in the example data set, ex.data
#' adx(x=ex.data[,1],g1=ex.data[,6],g2=ex.data[,7],d=.5)
#' 
#' #Standardized AD(x) for subpopulation 4 and subpopulation 5 in the example data set, ex.data
#' adx(x=ex.data[,1],g1=ex.data[,6],g2=ex.data[,7],d=.5,s=4.2)
#' @return A data frame of AD(x) indices, conditioned on the score scale
#' @return A plot of the AD(x) indices in reference to the difference that matters
#' @export

adx<- function(x,g1,g2,d,s,ymax,xlab,color){
  if(missing(s))
    s <- 1 
  
  ad <- (sqrt((g1-g2)^2))/s  
  dfdx <- data.frame("x"=x,"adx"=ad)
  
  if(missing(ymax))
    ymax <- max((dfdx[,2])+3)
  if(missing(xlab))
    xlab <- "Score Scale"
  if(missing(color))
    color <- "red"
 
  
  plot<- plot(type="b",x=dfdx[,1],y=dfdx[,2],ylim=c(0,ymax),xlab=xlab,ylab="Absolute Score Difference",col=color,pch=19)
  abline(d/s,0,col="black",lty=2)
    legend(x=min(dfdx[,1]),y=ymax,c("AD(x)","DTM"),col=c(color,"black"),pch=c(19,-1),lty=c(0,5))
  
  return(dfdx)
  return(plot)
  
  
}