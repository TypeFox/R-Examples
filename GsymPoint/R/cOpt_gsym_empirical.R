cOpt_gsym_empirical <-
function(t, parameters)
{
  cte = parameters[1]
  x0b = parameters[2:length(parameters)]
  ROCt = 1-cdf_empirical_dist(x0b, 1-t)
  res <- 1-ROCt-cte*t
  return(res)
}
