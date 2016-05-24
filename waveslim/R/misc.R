my.acf <- function(x)
{
  n <- length(x)
  x <- c(x, rep(0, n))
  Re(fft(Mod(fft(x))^2, inverse = TRUE)/2/n^2)[1:n]
}

my.ccf <- function(a, b) {
  n <- length(a)
  a <- c(a, rep(0, n))
  b <- c(b, rep(0, n))
  x <- Re(fft(fft(a) * Conj(fft(b)), inverse=TRUE))/2/n^2
  x[c((n+2):(2*n), 1:n)]
}
