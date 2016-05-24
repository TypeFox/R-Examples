# Modifed from pnormGC() in tigerstats package
# setting graph default value to TRUE, no area returned to console, round mean and sd

normGraph <- function(bound,region="below",mean=0,sd=1,graph=TRUE) {
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
    axis(1,at=places,labels=c("",as.character(round(places[2:3],3)),""))
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
    axis(1,at=places,labels=c("",as.character(round(places[2:3],3)),""))
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
  
  # return(area) we need only the graphs

}#end of pnormGC