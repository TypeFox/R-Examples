setwelch<-function(X, win = min(80, floor(length(X)/10)), inc = min(24,
    floor(length(X)/30)), coef = 64, wintaper=0.05  )
{
  #########  this is code borrowed from stft (e1071)

  ##########  wincoef = taper
  
  if(missing(wintaper)) { wintaper=0.05 }
  
  numcoef <- 2 * coef
  if (win > numcoef) {
    win <- numcoef
    cat("setwelch: window size adjusted to", win, ".\n")
  }

   

  wincoef=applytaper( rep(1,win), p=wintaper)

    
  numwin <- trunc((length(X) - win)/inc)
   
  z <- matrix(0, numwin + 1, numcoef)
    y <- z
    st <- 1
    for (i in 0:numwin) {
        z[i + 1, 1:win] <- X[st:(st + win - 1)] * wincoef
        y[i + 1, ] <- fft(z[i + 1, ])
        st <- st + inc
    }
    Y <- list(values = Mod(y[, 1:coef]), windowsize = win, increment = inc, wintaper=wintaper)
    return(Y)
}
