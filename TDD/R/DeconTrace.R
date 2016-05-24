DeconTrace = function(x, DPZ, fl = 0.05, fh = NaN, bitweight = NULL, dec = 1){
  # bitweight is optional A-to-D constant (when data are in counts) [V/count]
  dt = DPZ$dt

  # oversample if decimation is used
  if(dec != 1){
    n = length(x)
    x = Oversample(x, dec)
  }
  x = x - mean(x)

  # set up an "inverse response" iDPZ
  iZpg = DPZ$Zpg
#  iDPZ$poles = DPZ$zeros
#  iDPZ$np = DPZ$nz
#  iDPZ$zeros = DPZ$poles
#  iDPZ$nz = DPZ$np
#  iDPZ$Knorm = 1/DPZ$Knorm
#  iDPZ$Sense = 1/DPZ$Sense
  iZpg$pole = DPZ$Zpg$zero
  iZpg$zero = DPZ$Zpg$pole
  iZpg$gain = 1/DPZ$Zpg$gain
 
  # convert counts to volts, if necessary
  if(!is.null(bitweight)){
#    iDPZ$Sense = iDPZ$Sense * bitweight
      iZpg$gain = iZpg$gain * bitweight
  }
  
  # obtain coefficients
#  ba = PZ2Coef(iDPZ, dt)

  # deconvolve to velocity
#  v = filter(ba$b, ba$a, x)
  

  if(!is.na(fl) && (fl > 0 && fl < 1/(2*dt))){
    fZpg = as.Zpg(butter(2, fl*2*dt, 'high'))
    iZpg$pole = c(iZpg$pole, fZpg$pole)
    iZpg$zero = c(iZpg$zero, fZpg$zero)
    iZpg$gain = iZpg$gain * fZpg$gain
#    v = filter(ba$b, ba$a, v)
  }
  if(!is.na(fh) && (fh > 0 && fh < 1/(2*dt))){
#    ba = butter(2, fh*2*dt, 'low')
#    v = filter(ba$b, ba$a, v)
    fZpg = as.Zpg(butter(2, fh*2*dt, 'low'))
    iZpg$pole = c(iZpg$pole, fZpg$pole)
    iZpg$zero = c(iZpg$zero, fZpg$zero)
    iZpg$gain = iZpg$gain * fZpg$gain
  }

  # filter it
  v = filter(iZpg, x)
  # decimate if necessary
  if(dec != 1){
    v = v[1:n * dec]
  }
  return(v)
}
