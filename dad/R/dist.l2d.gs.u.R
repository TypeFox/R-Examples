dist.l2d.gs.u <-
function(x1, x2, check=FALSE)  
{
 # x1, x2  vectors containing the two samples
 if(check)
   {if(var(x1)<.Machine$double.eps | var(x2)<.Machine$double.eps) 
      {stop("At least one variance is zero") 
      }
   }
 return(sqrt(l2d.gs.u(x1, x1) + l2d.gs.u(x2, x2) - 2*l2d.gs.u(x1, x2)))
}

