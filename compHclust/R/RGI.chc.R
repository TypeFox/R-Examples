`RGI.chc` <-
function(x,xhat) {
  xbar <- rowMeans(x)
  return(rowSums((xhat-xbar)^2)/rowSums((x-xbar)^2))
}

