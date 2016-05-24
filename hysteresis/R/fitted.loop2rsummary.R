fitted.loop2rsummary <- function(object,...){
  g <- object
  return(data.frame("input"=g$pred.x,"output"=g$pred.y))
}
