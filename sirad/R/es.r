es <-
  function(Tmax,Tmin) {
    es <- (0.6108*exp(17.27*Tmax/(Tmax+237.3))+0.6108*exp(17.27*Tmin/(Tmin+237.3)))/2
    es
  }