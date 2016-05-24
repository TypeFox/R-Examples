shrinkageIntensity <-
function(hatThetaML = hatThetaML, n = n,
                               shrinkageTarget = shrinkageTarget){
  ans <- 1 - sum(hatThetaML^2)
  ans <- ans/((n - 1)*sum((shrinkageTarget - hatThetaML)^2))
  
  if(ans > 1) ans <- 1
  if(ans < 0) ans <- 0
  #check when n = c(0, 1)
  return(ans)
}
