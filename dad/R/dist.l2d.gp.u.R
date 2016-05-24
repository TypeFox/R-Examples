dist.l2d.gp.u <-
function(mean1, var1, mean2, var2, check=FALSE)  
{
  # mean1, var1 :   mean and variance of the first Gaussian density
  # mean2, var2 :   mean and variance of the second Gaussian density
  if(check)
  {if(abs(var1)<.Machine$double.eps | abs(var2)<.Machine$double.eps)
    {stop("At least one variance is zero")
    }
  }
  return(sqrt(l2d.gp.u(mean1, var1, mean1, var1) + l2d.gp.u(mean2, var2, mean2, var2) - 2*l2d.gp.u(mean1, var1, mean2, var2)))
}

