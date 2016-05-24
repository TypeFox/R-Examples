#' @title Graphical Calculator for the Empirical Rule

#' @description An app to facilitate visual understanding of Empirical Rule approximations
#' of probabilities, percentages.
#' 
#' @rdname EmpRuleGC
#' @usage EmpRuleGC(mean=0,sd=1,xlab="x")
#' @param mean Mean of the distribution
#' @param sd  Standard deviation of the distribution
#' @param xlab x-axis label  
#' @return Returns no value.  Used for the plotting side-effects.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @note  Uses \code{manipulate} in RStudio
#' @examples
#' \dontrun{
#' if(require(manipulate)) EmpRuleGC(mean=70,sd=3,xlab="Height (inches)")
#' }
EmpRuleGC <- function(mean=0,sd=1,xlab="x") {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  manipulate(
    range=picker("mean +- one SD"=1,
                 "mean +- two SDs"=2,
                 "mean +- three SDs"=3),
    shade=picker("between","outside"),
{curve(dnorm(x,mean=mean,sd=sd),from=mean-4*sd,to=mean+4*sd,xlab=xlab,ylab="density",
       main=paste("Empirical Rule:\nmean =",mean,", SD =",sd),axes=FALSE)
 DesiredLabels <- rep("",9)
 DesiredLabels[5] <- as.character(mean)
 DesiredLabels[5-range] <- as.character(mean-sd*range)
 DesiredLabels[5+range] <- as.character(mean+sd*range)
 axis(side=1,at=mean+sd*seq(-4,4,by=1),labels=DesiredLabels)
 
 if (shade=="between")  {
   x.coords <- seq(from=mean-sd*range,to=mean+sd*range,
                   length.out=500)
   y.coords <- dnorm(x.coords,mean=mean,sd=sd)
   x.coords <- c(mean-sd*range,x.coords,mean+sd*range)
   y.coords <- c(0,y.coords,0)
   polygon(x.coords,y.coords,col="lightblue")
   percs <- c("~68%","~95%","~99.7%")
   text(x=mean,y=0.45/(sqrt(2*pi)*sd),labels=percs[range],cex=2)
 } else {
   x.coords.left <- c(mean-4*sd,
                      seq(from=mean-4*sd,to=mean-sd*range,length.out=300),
                      mean-sd*range)
   x.coords.right <- c(mean+sd*range,
                       seq(from=mean+sd*range,to=mean+4*sd,length.out=300),
                       mean+4*sd)
   y.coords.left <- c(0,dnorm(seq(from=mean-4*sd,to=mean-sd*range,length.out=300),mean=mean,sd=sd),0)
   y.coords.right <- c(0,
                       dnorm(seq(from=mean+sd*range,to=mean+4*sd,length.out=300),mean=mean,sd=sd),
                       0)
   polygon(x.coords.left,y.coords.left,col="lightblue")
   polygon(x.coords.right,y.coords.right,col="lightblue")
   x.txt <- mean-sd*c(2.3,3,3.55)
   y.txt <- c(0.4,0.25,0.15)/(sqrt(2*pi)*sd)
   percs <- c("~16%","~2.5%","~0.15%")
   text(x.txt[range],y.txt[range],percs[range],cex=2)
   x.txt.right <- mean+sd*c(2.3,3,3.55)
   text(x.txt.right[range],y.txt[range],percs[range],cex=2)
                            
 }
    }
    )  #end manipulate
}  #end EmpRuleGC

if(getRversion() >= "2.15.1")  utils::globalVariables(c("shade"))
