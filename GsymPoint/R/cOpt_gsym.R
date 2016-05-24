cOpt_gsym <-
function(t,parameters,model)
{
  if (model == "norm")
  {
    # Under the binormal and bilognormal models:
    a = parameters[1]
    b = parameters[2]
    rho = parameters[3]
    val = pnorm((qnorm(1-rho*t)-a)/b)
    dif = val-t
  }
  res <- dif
  return(res)
}
