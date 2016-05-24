predicted.y <-
function(data, coeffs)
{
  
  num.coeffs = length(coeffs)
  
  pred = NULL
  
  for(i_data in 1:nrow(data))
  {
    y = coeffs[1] + sum(data[i_data,] * coeffs[-1])
    pred = c(pred , y)
  }
  
  return(as.numeric(pred))
}
