transfft <- function(sn, top=0.3){
  sn.fft = fft(sn)
  #plot(Re(sn.fft), type = "l")
  #qq <- round(length(sn)*.08)
  qq <- length(sn)*top # 3000
  #plot(Re(sn.fft), xlim = c(0, 1000), type = "l")
  sn.fft[qq:(length(sn)-qq)] = 0 + 0i
  sn.ifft = fft(sn.fft, inverse = TRUE)/length(sn.fft)
  return(Re(sn.ifft))
} 
