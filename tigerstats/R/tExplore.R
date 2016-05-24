#' @title Exploration of the t-Distributions

#' @description Plot the density curve of a t random variable at various degrees of freedom.  Compare
#' with the standard normal curve.
#' 
#' @rdname tExplore
#' @usage tExplore()
#' @return Used only for graphical side effects.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' if (require(manipulate)) tExplore()
#' }
tExplore <- function() {
  
  if (!("manipulate"  %in% installed.packages())) {
    return(cat(paste0("You must be on R Studio with package manipulate installed\n",
                      "in order to run this function.")))
  }
  
  manipulate(
    df=slider(1,100,step=1,initial=1,label="Degrees of Freedom"),
    ShowNorm=checkbox(FALSE,"Show Standard Normal Curve"),
{
  curve(dt(x,df=df),from=-4,to=4,ylab="density",n=1001,xlab="t",ylim=c(0,dnorm(0,0,1)),
        main=paste("t Density Curve, df = ",df))
  if (ShowNorm) curve(dnorm(x,0,1),from=-4,to=4,ylab="",n=1001,xlab="",col="red",lwd=2,
                    add=TRUE)
}    
    ) #end manipulate
} #end tExplore

if(getRversion() >= "2.15.1")  utils::globalVariables(c("ShowNorm"))
