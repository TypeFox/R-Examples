mcc <-
function(x, y, coeffs)
{
  y.hat = predicted.y(x, coeffs)
  mean.y.hat = mean(y.hat)
  res.y.hat  = y.hat - mean.y.hat
  
  mean.y = mean(y)
  res.y = y - mean.y
  
  numerator = sum(res.y * res.y.hat)
  
  denominator = sqrt(sum(res.y^2) * sum(res.y.hat^2))
  
  R = numerator / denominator
  
  return(R)
}
