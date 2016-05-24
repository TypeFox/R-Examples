l2d.gp <-
function(mean1,var1,mean2,var2,check=FALSE)
{
 if(check)
  {if(abs(det(var1))<.Machine$double.eps | abs(det(var2))<.Machine$double.eps)
    {stop("One of the sample variances is degenerate")
    }
  }
  p=length(mean1);
  d=mean1-mean2;
  vars=var1+var2;
  return(as.numeric((1/(2*pi)^(p/2))*(1/det(vars)^(1/2))*exp((-1/2)*t(d)%*%solve(vars)%*%d)))
}

