residuals.fittedloop <- function(object,...){
  g <- object
  resid.x <-g$pred.x-g$x
  resid.y <-g$pred.y-g$y
  resid.geometric <- sqrt(resid.x^2+resid.y^2)
  return(data.frame("input"=resid.x,"output"=resid.y,"geometric"=resid.geometric))
}
