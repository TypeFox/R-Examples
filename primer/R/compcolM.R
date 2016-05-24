`compcolM` <-
  function(t, y, params)
{
## This models S species with a competition-colonization tradeoff
  ## This requires 'params' is a list with named elements
  S <- params[["S"]]; D <- params[["D"]]
  with(params,
list( dpi.dt <- sapply(1:S, function(i) {
params[["ci"]][i] * y[i] * ( 1 - D - sum(y[1:i]) ) - 
   params[["m"]][i] * y[i] - sum( params[["ci"]][0:(i-1)] * y[0:(i-1)] * y[i] )  }
                       )
     )
       )
}
