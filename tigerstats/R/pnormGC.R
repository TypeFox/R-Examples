#' @title Graphical Calculator for Normal Curve Probabilities

#' @description Shades desired areas under a specified normal curve, returns numerical value of the area.
#' 
#' @rdname pnormGC
#' @usage pnormGC(bound,region="below",mean=0,sd=1,graph=FALSE)
#' @param bound A numerical vector of length 1 or 2, indicating boundary(ies) of shaded region on horizontal axis
#' @param region A character string.  Default is "below".  Possible values are "between" (when boundary consists of two numbers),
#' "below", "above", and "outside" (again when boundary consists of two numbers)
#' @param mean Mean of the distribution
#' @param sd  Standard deviation of the distribution 
#' @param graph Will produce graph of the probability
#' @return Numerical value of area under curve over region.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #This gives P(X < 75) for X normal with mean=70 and sd=4:
#' pnormGC(75,region="below",mean=70,sd=4)
#' 
#' #This gives P(X > 71) for X normal with mean=70 and sd=4:
#' pnormGC(71,region="above",mean=70,sd=4)
#' 
#' #This gives P(-1 < X < 1), for standard normal X:
#' pnormGC(c(-1,1),region="between")
#' 
#' #This gives P(X < 68 or X > 71), for X normal with mean =70 and sd=4:
#' pnormGC(c(68,71),region="outside",mean=70,sd=4)
pnormGC <- function(bound,region="below",mean=0,sd=1,graph=FALSE) {
  if (!is.numeric(bound)) stop("Specify one or two numerical boundaries")
  if (length(bound)==1 & !(region %in% c("below","above"))) stop("Specify region=\"below\" or
          region=\"above\"")
  if (length(bound)==2 & !(region %in% c("between","outside"))) stop("Specify region=\"between\" or
          region=\"outside\"")
  if (length(bound)>2) stop("Specify one or two numerical boundaries")

  if (length(bound)==2 & bound[1]>bound[2])  bound <- rev(bound)

  if (grepl("^be[lf]",region,perl=TRUE))  {
    area <- pnorm(bound,mean=mean,sd=sd)
    if (graph) {
    upper <- max(qnorm(.9999,mean=mean,sd=sd),bound+0.1*sd)
    lower <- min(qnorm(0.0001,mean=mean,sd=sd),bound-0.1*sd)
    curve(dnorm(x,mean=mean,sd=sd),from=lower,to=upper,ylab="density",axes=FALSE,n=50,
          main=paste("Normal Curve, mean = ",round(mean,2),", SD = ",round(sd,2),"\n Shaded Area = ",
                     round(area,4)))
    UnderShade(low=lower,high=bound,func=dnorm,mean=mean,sd=sd)
    axis(2)
    places <- c(lower,bound,mean,upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }

  if (grepl("^a[bf]",region,perl=TRUE))  {
    area <- pnorm(bound,mean=mean,sd=sd,lower.tail=FALSE)
    if (graph) {
    upper <- max(qnorm(.9999,mean=mean,sd=sd),bound+0.1*sd)
    lower <- min(qnorm(0.0001,mean=mean,sd=sd),bound-0.1*sd)
    curve(dnorm(x,mean=mean,sd=sd),from=lower,to=upper,ylab="density",axes=FALSE,n=50,
          main=paste("Normal Curve, mean = ",round(mean,2),", SD = ",round(sd,2),"\n Shaded Area = ",
                     round(area,4)))
    UnderShade(low=bound,high=upper,func=dnorm,mean=mean,sd=sd)
    axis(2)
    places <- c(lower,bound,mean,upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  if (grepl("^bet|^in",region,perl=TRUE))  {
    area <- pnorm(bound[2],mean=mean,sd=sd)-pnorm(bound[1],mean=mean,sd=sd)
    if (graph) {
    upper <- max(qnorm(.9999,mean=mean,sd=sd),bound+0.1*sd)
    lower <- min(qnorm(0.0001,mean=mean,sd=sd),bound-0.1*sd)
    curve(dnorm(x,mean=mean,sd=sd),from=lower,to=upper,ylab="density",axes=FALSE,n=50,
          main=paste("Normal Curve, mean = ",round(mean,2),", SD = ",round(sd,2),"\n Shaded Area = ",
                     round(area,4)))
    UnderShade(low=bound[1],high=bound[2],func=dnorm,mean=mean,sd=sd)
    axis(2)
    places <- c(lower,bound[1],bound[2],upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  if (grepl("^out",region,perl=TRUE))  {
    area <- pnorm(bound[1],mean=mean,sd=sd)+pnorm(bound[2],mean=mean,sd=sd,lower.tail=FALSE)
    if (graph) {
    upper <- max(qnorm(.9999,mean=mean,sd=sd),bound+0.1*sd)
    lower <- min(qnorm(0.0001,mean=mean,sd=sd),bound-0.1*sd)
    curve(dnorm(x,mean=mean,sd=sd),from=lower,to=upper,ylab="density",axes=FALSE,n=50,
          main=paste("Normal Curve, mean = ",round(mean,2),", SD = ",round(sd,2),"\n Shaded Area = ",
                     round(area,4)))
    UnderShade(low=lower,high=bound[1],func=dnorm,mean=mean,sd=sd)
    UnderShade(low=bound[2],high=upper,func=dnorm,mean=mean,sd=sd)
    axis(2)
    places <- c(lower,bound[1],bound[2],upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  return(area)

}#end of pnormGC
