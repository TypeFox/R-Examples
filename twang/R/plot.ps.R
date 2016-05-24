plot.ps <- function(x,plots="optimize",subset=NULL, color = TRUE,...)
{
   # Creates diag.plot plots and sends to current device
   # x:     ps object 
   # label: Label added to the plot titles

   # extract the propensity scores and weights from the ps object

pt1 <- diag.plot.color(x, plots, subset = subset, color = color, ...)

return(pt1)

}
