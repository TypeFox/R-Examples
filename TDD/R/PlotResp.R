PlotResp = function(PZ, DPZ, fmin = 0.01){
#  require(signal)
  N = 1/(DPZ$dt * fmin)
  f = (1:N - 1) * fmin
  trueresp = abs(PZ2Resp(PZ, f, FALSE))
#  coef = PZ2Coef(DPZ, DPZ$dt)
#  a = coef$a
#  b = coef$b
  impulse = 1:N == 1
  testresp = abs(fft(filter(DPZ$Zpg, impulse)))
  
  w = which(f > 0 & f < 0.5/DPZ$dt)
  plot(f[w], trueresp[w], type = 'l', log = 'xy')
  lines(f[w], testresp[w], col = 'red')
}
