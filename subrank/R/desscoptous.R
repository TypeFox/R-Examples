desscoptous <-
function(copest,normalize=FALSE)
{
  dim=length(copest$varnames)
  op <- par(mfrow=c(dim,dim))
  for (i in 1:dim)
  {
    for (j in 1:dim)
    {
      if (i==j)
      {
        plot.new()
        text(0.5,0.5,copest$varnames[i])
      }
      else
      {
        desscop(copest,copest$varnames[i],copest$varnames[j],normalize,FALSE)
      }
    }
  }
  par(op)
}

