`hilbert` <-
function(x)
  {
### calculate the hilbert transform of a signal
    
###    nn = klen/2;
###    numfreqs2 = 1+klen/2;
###    no2 = klen/2;
###    numfreqs2 = 1+klen/2;

    n = length(x)
    
    ff = fft(x)
    h=rep(0, n)
    
    if(n>0 & 2*floor(n/2)==n)
      {
### even and nonempty
        h[c(1, n/2+1)] = 1;
        h[2:(n/2)] = 2;
      }
    else
      {
        if(n>0)
          {
            ##  odd and nonempty
            h[1] = 1;
            h[2:((n+1)/2)] = 2;
          }
      }
    ht = fft(ff*h , inverse=TRUE)/length(ff);

    return(ht);

  }

