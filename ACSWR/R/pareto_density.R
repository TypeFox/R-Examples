pareto_density <-
function(x,scale,shape) {
  lpd <- ifelse(x<scale, -Inf, log(shape) + shape*log(scale) - (shape+1)*log(x))
  return(exp(lpd))
}
