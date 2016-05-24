MakeDPZ = function(PZ, dt, fmin = 1/360, niter = 50000, ...){
  # This function tests both the bilinear method and the finite difference approximation method, and returns the results of the better approximation.

    # finite difference:
  N = 2^ceiling(log(1/(fmin * dt), 2)) # to make sure the fft is fast
  imp = 1:N == 1
  Xf = MatchCoefDPZ(PZ, dt, N, niter = niter, ...)$DPZ

  # bilinear:
  Xb = PZ
  Xb$dt = dt
  Xb$Zpg = bilinear(Sz = PZ$zeros, Sp = PZ$poles, Sg = PZ$Sense * PZ$Knorm, T = dt)

  f = (1:N - 1)/(N * dt)
  trueresp = abs(PZ2Resp(PZ, f, PLOT = TRUE))
  fresp = abs(fft(ConvolveTrace(imp, Xf)))
  bresp = abs(fft(ConvolveTrace(imp, Xb)))

  w = which(f < 0.5/dt)
  lines(f[w], fresp[w], col = 'red')
  lines(f[w], bresp[w], col = 'blue')

  legend('topleft', lty = 1, col = c(1, 2, 4), legend = c('True', 'FD', 'Bilinear'))
  ff = f[which(abs((fresp - trueresp)/trueresp) > 0.01 & f > f[5] & f < 0.5/dt)[1]]
  if(is.na(ff)){
      ff = 0.5/dt
  }
  fb = f[which(abs((bresp - trueresp)/trueresp) > 0.01 & f > f[5] & f < 0.5/dt)[1]]
  if(is.na(fb)){
      fb = 0.5/dt
  }
  print(paste('finite-difference:', ff, ';', 'bilinear transform:', fb))
  if(ff > fb){
      Xf$fmax = ff
      Xf$method = 'FD'
      return(Xf)
  }else{
      Xb$fmax = fb
      Xb$method = 'bilinear'
      return(Xb)
  }
}
