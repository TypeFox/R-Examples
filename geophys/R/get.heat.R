get.heat <-
function(x, T0, k, t)
  {
	tim = t
        t1 = erf(x/(2*sqrt(k*tim)))
        T = T0*(0.5+0.5*c(t1))
    return(T)
  }

