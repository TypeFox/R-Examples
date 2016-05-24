residuals.ellipsefit <- function(object,...){
  g <- object
  resid.x <-g$pred.x-g$x
  resid.y <-g$pred.y-g$y
  resid.geometric <- sqrt(resid.x^2+resid.y^2)
  resid.algebraic <- g$residuals
  return(data.frame("input"=resid.x,"output"=resid.y,"geometric"=resid.geometric,"algebraic"=resid.algebraic))
}
