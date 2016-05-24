jitter.lab <-
function (x,w)
{

  ox = order(x)
  x = x[ox]
  
  L=length(x)
  y=rep(1,L)
  if(length(w)==1){w=rep(w,L)}
  for(i in 1:L)
    {
      b=which(x<x[i] & (w+x)>=x[i])
      y[i]=min(which(!(1:L %in% y[b])))
    }
  y = y[order(ox)]
  return(y-1)

}

