`RESCALE` <-
function(x, nx1, nx2, minx, maxx)
{
  #    rescale a vector 
  nx = nx1+(nx2-nx1)*(x-minx)/(maxx-minx)
  return(nx)
}

