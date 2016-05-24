`scal2freqs` <-
function(octs, dt, plot=FALSE)
{
  ####   make a vector of  freqencies from the scales
  if(missing(plot) ) { plot=FALSE }
  mm  = Mmorlet(-8, 8, 256)
   ##  m2 = morlet(256, 128,  256/16, w0=5)
  psi = mm$morl
  psiFT = fft(psi); 
  sp = (abs(psiFT)); 
  indmax = which.max(sp)
  vmax = sp[indmax]
  TD = max(mm$xval)-min(mm$xval);
  per = TD/(indmax-1); 
  freq = 1/per;

  if(plot==TRUE)
    {
      psiFT[sp<vmax] = 0;
      
      recfreq = fft(psiFT, inverse = TRUE);
      plot(mm$xval, mm$morl, type='l')
      
      lines(mm$xval, 0.75*max(abs(psi))*Re(recfreq)/max(abs(recfreq)), col=2)
    }

#####
  freqs = freq / (octs * dt)


   return(freqs)

}

