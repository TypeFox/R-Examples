FourierRand <-
function(x){
  n <- length(x)
  z <- fft(x)

  if(n%%2 == 0){
    ph <- 2*pi*runif(n/2-1)
    ph <- c(0, ph, 0, -rev(ph))}

  if(n%%2 != 0){
    ph <- 2*pi*runif((n-1)/2)
    ph <- c(0, ph, -rev(ph))}

  ph <- complex(imaginary = ph)
  z <- z * exp(ph)
  x.sur <- Re(fft(z, inverse = TRUE)/n)
  
  return(invisible(x.sur))
}

## Code: Tian, H. and Cazelles, B., \code{WaveletCo}
