#' Root Mean Square Difference
#' @param x a column vector of scores on which the rmsd is conditioned
#' @param o a column vector of equated scores based on the overall population (aligned with elements in x)
#' @param g column vectors of equated scores based on various subpopulations (aligned with elements in x)
#' @param w a row vector of weights for subpopulations 1 thru n (length = number of groups)
#' @param d a scalar of the difference that matters
#' @param s a scalar representing the standard deviation of x for any (sub)population of interest (e.g., synthetic population) (default is 1, which leads to calculation of the unstandardized rmsd)
#' @param ymax A maximum value for the y axis of the RMSD graph (default = 3 + the maximum RMSD value)
#' @param xlab A label for the x axis of the RMSD graph (default = Score Scale)
#' @param color color of RMSD line and points on plot (default = red)
#' @author Anne Corinne Huggins-Manley 
#' @references 
#' \itemize{
#' \item{Dorans, N.J., & Holland, P.W. (2000). Population invariance and the equitability of tests: Theory and the linear case. Journal of Educational Measurement, 37, 281-306.}
#' }
#' @seealso \code{\link{remsd}}
#' @examples
#' #Unstandardized RMSD for subpopulations 1 and 2 in the example data set, ex.data, assuming equal weights for the subpopulations
#' rmsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4]),w=c(.5,.5),d=.5)
#' 
#' #Unstandardized RMSD for all five subpopulations in the example data set, ex.data
#' rmsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),w=c(.1,.2,.4,.2,.1),d=.5)
#' 
#' #Unstandardized RMSD for all five subpopulations in the example data set, ex.data, with adjustments to the maximum y-axis on the plot, a new xlabel, and points/line in blue
#' rmsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),w=c(.1,.2,.4,.2,.1),d=.5,ymax=3,xlab="Old Form",color="blue")
#' 
#' #Standardized RMSD for all five subpopulations in the example data set, ex.data
#' rmsd(x=ex.data[,1],o=ex.data[,2],g=c(ex.data[,3],ex.data[,4],ex.data[,5],ex.data[,6],ex.data[,7]),w=c(.1,.2,.4,.2,.1),d=.5,s=4.2)
#' @note The equally weighted version of this index (Kolen & Brennan, 2004) can be obtained by inputting a w vector consisting of identical elements that sum to 1. See example 1 above.
#' @return A data frame of RMSD indices, conditioned on the score scale
#' @return A plot of the RMSD in reference to the difference that matters
#' @export

rmsd<- function(x,o,g,w,d,s,ymax,xlab,color){
  
  if(missing(s))
    s <- 1
  
  tm<- c()  
  
  for (i in 1:length(w)){
    t<- assign(x=paste("t",i,sep=""),value=c(rep(w[i],length(x))))
    
    tm <- c(tm,t)  
  }
  
  
  rmsdX<-(tm*(o-g)^2)
  rmsdY<-rowSums(matrix(rmsdX,nrow=length(x))) 
  rmsdZ<-(sqrt(rmsdY))/s
  dfrmsd <- data.frame("x"=x,"rmsd"=rmsdZ)
  
  if(missing(ymax))
    ymax <- max((dfrmsd[,2])+3)
  if(missing(xlab))
    xlab <- "Score Scale"
  if(missing(color))
    color <- "red"
  
  plot<- plot(type="b",x=dfrmsd[,1],y=dfrmsd[,2],ylim=c(0,ymax),xlab=xlab,ylab="RMSD",col=color,pch=19)
  abline(d/s,0,col="black",lty=2)
  legend(x=min(dfrmsd[,1]),y=ymax,c("RMSD","DTM"),col=c(color,"black"),pch=c(19,-1),lty=c(0,5))
  
  return(dfrmsd)
  return(plot)
  
}