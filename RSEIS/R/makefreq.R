`makefreq` <-
function(n,dt)
  {
##% Make frequency vector for Fourier Transforms
##% USAGE: f=makefreq(n,dt);
##%  construct the frequency vector that is consistent with the result of an FFT.
##%  dt is the time-domain sample interval and n is the number of points.
##%
##%  See also FT and IFT.

    if( (n%%2) == 1 )
      {
        f=c(seq(0,(n-1)/2),   seq(-(n-1)/2, -1))/ (n*dt) 
                                        #   frequency vector for odd n
      }
    else
      {
        f=c( seq(0,n/2)   ,     seq(-n/2+1,-1))/(n*dt);   #  frequency vector for even n
      }
    return(f)
  }

