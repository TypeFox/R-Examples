idcost <-
function( Y) { 

  d1 = dim(Y)
  n = d1[1]
  m  = d1[2]
  
  X=matrix(0,n,m)
  ci = complex(real = 0, imaginary = (0:(n-1)) )
  
  
###  Precompute weights
  ww = sqrt(2*n) *   exp(ci*pi/(2*n))

### loop over columns
  for (  i in 1 : m )
    {
      if( (n %% 2)  ==1)
        {
###  odd case
          if(i==1)
            {
              ww[1] = ww[1] * sqrt(2)
            }
          
          yy = rep(0,2*n)
          yy[1:n] = ww * Y[,i]
          yy[(n+2):(2*n)] = complex( imaginary=-1) *  ww[2:n]   * Y[seq(from=n, by=-1, to=2) ,i]

          y = fft(yy, inverse=TRUE)/length(yy)
        }
      else
        {
###  even case
          if(i==1){ ww[1] = ww[1]/sqrt(2)   }
          
          yy = ww * Y[,i]
          yt = fft(yy, inverse=TRUE)/length(yy)
          y = matrix(0,n,1)
          y[ seq(from=1, by=2, to=n) ] = yt[1:(n/2) ]
          y[seq(from=2, by=2, to=n) ]  = yt[seq(from=n, by=-1, to=(n/2+1) )]
        } 
      
###  Extract inverse DCT getting rid of the imaginary roundoff error
      X[,i] = Re(y[1:n])
      
    }

  return(X)

}
