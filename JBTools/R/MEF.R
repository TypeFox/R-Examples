MEF <-function(
  ##title<< Compute the modelling efficiency (MEF)
  prediction    ##<< numeric: array 1
  , observation ##<< numeric: array 2
)
  ##description<< Calculates the modelling efficiency (MEF) of two arrays 
{
  if (length(prediction) != length(observation))
    stop('Vectors need to of the same length!')
  
  valid <- !is.na(prediction) &  !is.na(observation)
  prediction <-  prediction[valid]
  observation  <-  observation[valid]

  
  if (length(prediction) == 0 | length(observation) == 0) {
    value_out = NaN
  } else {
    value_out=  (1 - (sum(( prediction  - observation) ^ 2)) / 
          (sum((observation - mean(observation)) ^ 2)))
  }     
  ##value<< MEF of the two vectors
  return(value_out)
}
