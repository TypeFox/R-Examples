ConvolveTrace = function(x, DPZ, dec = 1){
  if(dec != 1){
    n = length(x)
    x = Oversample(x, dec)
  }
  dt = DPZ$dt
  # obtain coefficients
  ba = PZ2Coef(DPZ, dt)

  # convolve to voltage
  volt = filter(DPZ$Zpg, x)

  # decimate if necessary
  if(dec != 1){
    volt = volt[1:n * dec]
  }
  return(volt)
}
