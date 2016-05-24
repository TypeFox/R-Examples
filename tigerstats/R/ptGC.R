#' @title Graphical Calculator for t-Curve Probabilities

#' @description Shades desired areas under a specified t-curve, returns numerical value of the area.
#' 
#' @rdname ptGC
#' @usage ptGC(bound,region="between",df=1,graph=FALSE)
#' @param bound A numerical vector of length 1 or 2, indicating boundary(ies) of shaded region on horizontal axis
#' @param region A character string.  Possible values are "between" (when boundary consists of two numbers),
#' "below", "above", and "outside" (again when boundary consists of two numbers)
#' @param df degrees of freedom of the distribution
#' @param graph produce graph?
#' @return Numerical value of area under curve over region.  Also plots the t-curve with the shaded area.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #This gives P(-2 < t < 2) for a t-random variable with 1 degree of freedom:
#' ptGC(c(-2,2),region="between",df=1)
#' 
#' #This gives P(t < -1) for a t-random variable with 5 degrees of freedom:
#' ptGC(-1,region="below",df=5)
#' 
#' #This gives P( t < -2 OR t >2), for a t-random variable with 5 degrees of freedom:
#' ptGC(c(-2,2),region="outside",df=5)
ptGC <- function(bound,region="between",df=1,graph=FALSE) {
  if (!is.numeric(bound)) stop("Specify one or two numerical boundaries")
  if (length(bound)==1 & !(region %in% c("below","above"))) stop("Specify region=\"below\" or
                                                                 region=\"above\"")
  if (length(bound)==2 & !(region %in% c("between","outside"))) stop("Specify region=\"between\" or
                                                                     region=\"outside\"")
  if (length(bound)>2) stop("Specify one or two numerical boundaries")
  
  if (length(bound)==2 & bound[1]>bound[2])  bound <- rev(bound)
  
  if (region=="below")  {
    area <- pt(bound,df=df)
    if (graph) {
    upper <- max(4,bound+1)
    lower <- min(-4,bound-1)
    curve(dt(x,df=df),from=lower,to=upper,ylab="density",axes=FALSE,n=50,ylim=c(0,1.1*dt(0,df=df)),
          main=paste("t-curve, df = ",df,"\n Shaded Area = ",round(area,4)))
    UnderShade(low=lower,high=bound,func=dt,df=df)
    axis(2)
    places <- c(lower,bound,0,upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  if (region=="above")  {
    area <- 1-pt(bound,df=df)
    if (graph) {
    upper <- max(4,bound+1)
    lower <- min(-4,bound-1)
    curve(dt(x,df=df),from=lower,to=upper,ylab="density",axes=FALSE,n=50,ylim=c(0,1.1*dt(0,df=df)),
          main=paste("t-curve, df = ",df,"\n Shaded Area = ",round(area,4)))
    UnderShade(low=bound,high=upper,func=dt,df=df)
    axis(2)
    places <- c(lower,bound,0,upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  if (region=="between")  {
    area <- pt(bound[2],df=df)-pt(bound[1],df=df)
    if (graph) {
    upper <- max(4,bound+1)
    lower <- min(-4,bound-1)
    curve(dt(x,df=df),from=lower,to=upper,ylab="density",axes=FALSE,n=50,ylim=c(0,1.1*dt(0,df=df)),
          main=paste("t-curve, df = ",df,"\n Shaded Area = ",round(area,4)))
    UnderShade(low=bound[1],high=bound[2],func=dt,df=df)
    axis(2)
    places <- c(lower,bound[1],bound[2],upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  if (region=="outside")  {
    area <- pt(bound[1],df=df)+pt(bound[2],df=df,lower.tail=FALSE)
    if (graph) {
    upper <- max(4,bound+1)
    lower <- min(-4,bound-1)
    curve(dt(x,df=df),from=lower,to=upper,ylab="density",axes=FALSE,n=50,ylim=c(0,1.1*dt(0,df=df)),
          main=paste("t-curve, df = ",df,"\n Shaded Area = ",round(area,4)))
    UnderShade(low=lower,high=bound[1],func=dt,df=df)
    UnderShade(low=bound[2],high=upper,func=dt,df=df)
    axis(2)
    places <- c(lower,bound[1],bound[2],upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  return(area)
  
}#end of ptGC
