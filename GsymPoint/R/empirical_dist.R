empirical_dist <-
function(x,t)
{
  x_logico = (x<=t)
  y = sum(x_logico)/length(x)
  res <- y
  return(res)
}
