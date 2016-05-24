deltaVP <-
  function(Tmax,Tmin) {
    Tmean <- (Tmax+Tmin)/2
    deltaVP <- (4098*(0.6108*exp(17.27*Tmean/(Tmean+237.3))))/(Tmean+237.3)^2
    deltaVP
  }