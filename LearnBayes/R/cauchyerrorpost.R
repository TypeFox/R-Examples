cauchyerrorpost=function(theta, data)
{
logf=function(data,theta)
  log(dt((data-theta[1])/exp(theta[2]),df=1)/exp(theta[2]))

return(sum(logf(data,theta)))
}
