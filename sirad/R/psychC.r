psychC <-
  function(Tmax,Tmin,z) {
    p <- 101.3*((293-0.0065*z)/293)^5.26
    psychC <- (1.013e-03*p)/(0.622*(2.501-0.002361*((Tmax+Tmin)/2)))
    psychC
  }