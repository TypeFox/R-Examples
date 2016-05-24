#' @title Graphical Calculator for Normal Curve Percentiles

#' @description When you know a certain area under a normal denisity durve, this function returns
#' the x-axis values of the boundary of that area.
#' 
#' @rdname qnormGC
#' @usage qnormGC(area,region="below",mean=0,sd=1,graph=FALSE)
#' @param area The known percentile
#' @param region A character string.  Default is "below".  Other possible values are "between" 
#' (when known area is symmetric around the mean two numbers),
#' "below", "above", and "outside" (when knonw area is outside a region symmetric around the mean)
#' @param mean Mean of the distribution
#' @param sd  Standard deviation of the distribution 
#' @param graph Will produce graph of the area
#' @return Numerical value of the percentile, and a vector when there are two bounds.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #80th percentile of a normal distribution with mean=70 and sd=4:
#' qnormGC(0.80,region="below",mean=70,sd=4)
#' 
#' #Return value x so that P(X > x) = 0.10 (same as the 90th percentile)
#' qnormGC(0.10,region="above",mean=70,sd=4)
#' 
#' #This gives the multiplier for 95%-confidence intervals based on the z-statistic
#' qnormGC(0.95,region="between")
#' 
#' #This gives critical values for a two-sided z-test with alpha = 0.01:
#' qnormGC(0.01,region="outside")
qnormGC <- function(area,region="below",mean=0,sd=1,graph=FALSE) {
  if (!is.numeric(area)) stop("Specify an area between 0 and 1")
  if (area >= 1) stop("Express area as a number between 0 and 1")
  if (length(area)==1 & !(region %in% c("below","above","between","outside"))) stop('Specify region="below" or
          "above" or "between" or "outside"')
  if (length(area)>1) stop("Specify just one area, please!")


  if (region=="below")  {
    quant <- qnorm(area,mean=mean,sd=sd)
    if (graph) {
    upper <- max(qnorm(.9999,mean=mean,sd=sd),quant+0.1*sd)
    lower <- min(qnorm(0.0001,mean=mean,sd=sd),quant-0.1*sd)
    curve(dnorm(x,mean=mean,sd=sd),from=lower,to=upper,ylab="density",axes=FALSE,n=50,
          main=paste0("Normal Curve, mean = ",round(mean,2),", SD = ",round(sd,2),"\n Percentile = ",
                     round(100*area,4),"%"))
    UnderShade(low=lower,high=quant,func=dnorm,mean=mean,sd=sd)
    axis(2)
    places <- c(lower,round(quant,2),upper)
    axis(1,at=places,labels=c("",as.character(places[2]),""))
    }
  }

  if (region=="above")  {
    quant <- qnorm(area,mean=mean,sd=sd,lower.tail=FALSE)
    if (graph) {
    upper <- max(qnorm(.9999,mean=mean,sd=sd),quant+0.1*sd)
    lower <- min(qnorm(0.0001,mean=mean,sd=sd),quant-0.1*sd)
    curve(dnorm(x,mean=mean,sd=sd),from=lower,to=upper,ylab="density",axes=FALSE,n=50,
          main=paste0("Normal Curve, mean = ",round(mean,2),", SD = ",round(sd,2),"\n Percent Above = ",
                     round(100*area,4),"%"))
    UnderShade(low=quant,high=upper,func=dnorm,mean=mean,sd=sd)
    axis(2)
    places <- c(lower,round(quant,2),upper)
    axis(1,at=places,labels=c("",as.character(places[2]),""))
    }
  }
  
  if (region=="between")  {
    areaBelow <- area + 0.5*(1-area)
    quantUpper <- qnorm(areaBelow,mean=mean,sd=sd)
    quantLower <- qnorm(1-areaBelow,mean=mean,sd=sd)
    quant <- c(quantLower,quantUpper)
    if (graph) {
    upper <- max(qnorm(.9999,mean=mean,sd=sd),quantUpper+0.1*sd)
    lower <- min(qnorm(0.0001,mean=mean,sd=sd),quantLower-0.1*sd)
    curve(dnorm(x,mean=mean,sd=sd),from=lower,to=upper,ylab="density",axes=FALSE,n=50,
          main=paste("Normal Curve, mean = ",round(mean,2),", SD = ",round(sd,2),"\n Percent Between = ",
                     round(100*area,4),"%"))
    UnderShade(low=quantLower,high=quantUpper,func=dnorm,mean=mean,sd=sd)
    axis(2)
    places <- c(lower,round(quantLower,2),round(quantUpper,2),upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  if (region=="outside")  {
    quantUpper <- qnorm(0.5*area,mean=mean,sd=sd,lower.tail=FALSE)
    quantLower <- qnorm(0.5*area,mean=mean,sd=sd)
    quant <- c(quantLower,quantUpper)
    if (graph) {
      upper <- max(qnorm(.9999,mean=mean,sd=sd),quantUpper+0.1*sd)
      lower <- min(qnorm(0.0001,mean=mean,sd=sd),quantLower-0.1*sd)
      curve(dnorm(x,mean=mean,sd=sd),from=lower,to=upper,ylab="density",axes=FALSE,n=50,
            main=paste("Normal Curve, mean = ",round(mean,2),", SD = ",round(sd,2),"\n Percent Outside = ",
                       round(100*area,4),"%"))
      UnderShade(low=quantLower,high=quantUpper,func=dnorm,mean=mean,sd=sd)
      axis(2)
      places <- c(lower,round(quantLower,2),round(quantUpper,2),upper)
      axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  return(quant)

}#end of qnormGC
