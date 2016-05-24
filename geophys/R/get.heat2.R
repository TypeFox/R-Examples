get.heat2 <-
function(x, T0, k, t)
  {
	tim = t
        t1 = erf(x/(2*sqrt(k*tim)))
        T = T0*(c(t1))
    return(T)
  }

