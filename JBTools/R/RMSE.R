RMSE <-function(
##title<< Compute the residual mean square error (RMSE)
   prediction    ##<< numeric: array 1
   , observation ##<< numeric: array 2
)
  ##description<< Calculates the residual mean square error (RMSE) of two arrays
{
  if (length(prediction) != length(observation))
    stop('Vectors need to of the same length!')
  if (sum((!is.na(prediction) & !is.na( observation))) == 0) {
    value_out = NaN
  } else {
    value_out=  sqrt(    sum((prediction- observation)^2) / length(prediction))
  }
  ##value<< RMSE of the two vectors
  return(value_out)
}
