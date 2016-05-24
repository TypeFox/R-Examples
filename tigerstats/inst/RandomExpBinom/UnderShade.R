#' @title Shade Under Density Curves

#' @description Utility function for ptGC, pnormGC, pchisqGC, possibly others
#' @keywords internal
#' @rdname UnderShade
#' @usage UnderShade(low,high,func,...)
#' @param low lower bound
#' @param high upper bound
#' @param func density function
#' @param \ldots other arguments passed (to modify func)
#' @return graphical side effect only
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
UnderShade <- function(low,high,func,...) { #Utility
  x.coords <- c(low,seq(low,high,length.out=301),high)
  y.coords <- c(0,func(seq(low,high,length.out=301),...),0)
  polygon(x.coords,y.coords,col="lightblue",cex=2)
}