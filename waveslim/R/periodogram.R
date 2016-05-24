per <- function(z) {
  n <- length(z)
  (Mod(fft(z))**2/(2*pi*n)) [1:(n %/% 2 + 1)]
}
