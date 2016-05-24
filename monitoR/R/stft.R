# Copied from seewave, with only modification being addition of fftw:: at all planFFT and FFT calls
# Modified: 2015 APR 2

stft <-
function (wave, f, wl, zp, step, wn, norm = FALSE, fftw = FALSE) 
{
   zpl <- zp%/%2
   if (zpl == 0) {
      W <- ftwindow(wl = wl, wn = wn)
      if (fftw) {
         p <- fftw::planFFT(wl)
         z1 <- apply(as.matrix(step), 1, function(x) Mod(fftw::FFT(wave[x:(wl + x - 1)] * W, plan = p)))
      }
      else {
         z1 <- apply(as.matrix(step), 1, function(x) Mod(fft(wave[x:(wl + x - 1), ] * W)))
      }
   }
   else {
      W <- ftwindow(wl = wl + zp, wn = wn)
      if (fftw) {
         p <- fftw::planFFT(wl + zp)
         z1 <- apply(as.matrix(step), 1, function(x) Mod(fftw::FFT(c(1:zpl, 
            wave[x:(wl + x - 1), ], 1:zpl) * W, plan = p)))
      }
      else {
         z1 <- apply(as.matrix(step), 1, function(x) Mod(fft(c(1:zpl, 
            wave[x:(wl + x - 1), ], 1:zpl) * W)))
      }
   }
   z2 <- z1[1:((wl + zp)/2), , drop = FALSE]
   if (norm) {
      z <- matrix(numeric(length(z2)))
      dim(z) <- dim(z2)
      for (i in 1:ncol(z)) {
         z[, i] <- z2[, i]/max(z2[, i])
      }
   }
   else {
      z <- z2/max(z2)
   }
   return(z)
}
