
StepSDE <- function (drift, diffusion, dt = 0.01) 
{
  sdt = sqrt(dt)
  return(function(x0, t0, deltat, ...) {
    x = x0
    t = t0
    termt = t0 + deltat
    v = length(x)
    repeat {
      dw = rnorm(v,0,sdt)
          x = x + drift(x, t, ...)*dt + as.vector(diffusion(x, t,...)%*%dw)
          t = t + dt
          if (t > termt)
            return(x)
        }
  })
}


# eof

