dist.l2d.gs <-
function(x1, x2, check=FALSE)  
{
  return(sqrt(l2d.gs(x1, x1, check = check) + l2d.gs(x2, x2, check = check) - 2*l2d.gs(x1, x2, check = check)))
}

