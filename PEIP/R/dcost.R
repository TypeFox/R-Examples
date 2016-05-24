dcost <-
function( X) {  

  if(is.matrix(X))
    {
      d1 = dim(X)
      n = d1[1]
      m = d1[2]
    }
  else
    {
      n = length(X)
      m  = 1
      X = matrix(X, n, m)
      
    }
  
  
  
  Y=matrix(0,n, m)
  
###  Precompute weights
  ci = complex(real = 0, imaginary = -(0:(n-1)) )
  
  
  ww = ( (exp(ci*pi/(2*n))/sqrt(2*n) )) 
  
  
  ww[1] = ww[1] / sqrt(2)
  
### loop over columns
  for (  i in 1:m )
    {
      if( (n %% 2) ==1)
        {
###  odd case
          y = rep(0,2*n)
          y[1:n] = X[ ,i]
          y[(n+1):(2*n) ] = X[seq(from=n, by=-1, to=1) ,i]
          ff = fft(y)  
          ff = ff[1:n]
        }
      else
        {
###  even case
          y = c(X[seq(from=1, by=2, to=n) ,i],  X[ seq(from=n, by=-2, to=2) ,i] )
          ff = fft(y)
          if(i==1){ ww = 2*ww   } 
        }
      
      
      Y[,i] = Re(ww * ff)
    } 
  
  return(Y)
  
}
