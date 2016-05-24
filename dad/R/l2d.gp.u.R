l2d.gp.u <-
function(mean1,var1,mean2,var2,check=FALSE)
{
  if(check)
  {if(abs(var1)<.Machine$double.eps | abs(var2)<.Machine$double.eps)
    {stop("At least one variance is zero")
    }
  }
  d=mean1-mean2;
  vars=var1+var2;
  return((1/sqrt(2*pi))*(1/sqrt(vars))*exp(-(1/2)*(d^2)/vars))
}

