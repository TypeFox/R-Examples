cred_int <-
function(x,y,xbins,conf.int=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95))
{
  if ( length(x)!=length(y) ) stop("x and y don't have the same length")
  
  bin <- hexbin(x,y, xbins=xbins)
  t <- rev(unique(sort(bin@count)))
  colorcut <- rep(NA,length(conf.int))
  s <- 0
  i <- 0
  for (j in 1:length(colorcut))
  {
    while ( s < sum(bin@count)*conf.int[j] )
    {
      i <- i+1
      s <- s + sum( bin@count[ bin@count == t[i] ] )
    }
    colorcut[j] <- (t[i]-0.5)/max(t)
  }
  
  list(t=t,i=i,s=s/sum(bin@count),colorcut=c(0,rev(colorcut),1),count=bin@count)
}