#' Comparative Plot of Two Distributions
#' @param cases name of case group data (vector)
#' 
#' @param control name of control group data (vector)
#' 
#' @return beanplot and Q-Q plot of comparing distribution pairs.
#' @author Yusuke Matsui & Teppei Shimamura
#' 
#' 
#' @examples
#' library(D3M)
#' cases <- rbeta(50,1,5)
#' control <- rbeta(50,2,5)
#' d3m.plot(cases,control)
#' 
#' @export
#' 

d3m.plot <- function(cases,control){

  #library(graphics)
  #library(beanplot)
  #library(stats)
  
  graphics::par(mar=c(2,2,2,1))
  
  graphics::par(mfrow = c(1,2))
  
  cases <- as.numeric(cases)
  
  control <- as.numeric(control)
  
  beanplot::beanplot(cases,control,side="both",col = list("black",c("grey","white")),ylim=c(-.1,1),ll=0.05,names = c("cases","control"),horizontal = T,log = "",main="Density of beta value")
  
  graphics::grid()
  
  stats::qqplot(cases,control,xlim=c(0,1),ylim=c(0,1),main="Q-Q plot")
  
  graphics::abline(a = 0,b = 1,col="red")
}
