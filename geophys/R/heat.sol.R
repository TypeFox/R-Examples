heat.sol <-
function(x, T0, k, t)
  {

    for(tim in t)
      {
        t1 = erf(x/(2*sqrt(k*tim)))
        t2 = rev(-t1)
        T = T0*(0.5+0.5*c(t2,t1))
        
        plot(c(rev(-x),x),T, xlab="Distance, m", ylab="Temp, K")
        abline(v=0, lty=2)
        atim = tim/(24*3600)
        
        title(paste(sep=' ', "time=",atim, " days") )
        locator()
      }
    return(T)
  }

