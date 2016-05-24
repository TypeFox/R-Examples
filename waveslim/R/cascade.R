wavelet.filter <- function(wf.name, filter.seq = "L", n = 512)
{
  cascade <- function(f, x, j)
    {
      L <- length(f)
      N <- length(x)
      M <- (L - 1) * 2^j
      M1 <- M - L + 2
      M2 <- 2 * M - L + 2
      if(N > M1)
        stop("x is too long\n")
      else x <- c(x, rep(0, M1 - N))
      xj <- c(rep(0, M), x, rep(0, M))
      yj <- rep(0, M2)
      for(i in 1:L)
        yj <- yj + f[L - i + 1] * xj[1:M2 + (i - 1) * 2^j]
      yj
    }
  if(is.character(wf.name))
    wf <- wave.filter(wf.name)
  else
    wf <- wf.name
  J <- nchar(filter.seq)
  key <- rev(substring(filter.seq, 1:J, 1:J))
  f <- 1
  fl <- wf$lpf
  fh <- wf$hpf
  for(k in 1:J) {
    if(key[k] == "H")
      f <- cascade(fh, f, k - 1)
    else if(key[k] == "L")
      f <- cascade(fl, f, k - 1)
    else stop("Invalid filter.seq\n")
  }
  f
}

squared.gain <- function(wf.name, filter.seq = "L", n = 512)
{
  cascade <- function(f, x, j)
    {
      L <- length(f)
      N <- length(x)
      M <- (L - 1) * 2^j
      M1 <- M - L + 2
      M2 <- 2 * M - L + 2
      if(N > M1)
        stop("x is too long\n")
      else x <- c(x, rep(0, M1 - N))
      xj <- c(rep(0, M), x, rep(0, M))
      yj <- rep(0, M2)
      for(i in 1:L)
        yj <- yj + f[L - i + 1] * xj[1:M2 + (i - 1) * 2^j]
      yj
    }
  if(is.character(wf.name))
    wf <- wave.filter(wf.name)
  else 
    wf <- wf.name
  J <- nchar(filter.seq)
  key <- rev(substring(filter.seq, 1:J, 1:J))
  f <- 1
  fl <- wf$lpf
  fh <- wf$hpf
  for(k in 1:J) {
    if(key[k] == "H")
      f <- cascade(fh, f, k - 1)
    else if(key[k] == "L")
      f <- cascade(fl, f, k - 1)
    else stop("Invalid filter.seq\n")
  }
  Mod(fft(c(f, rep(0, n - length(f))))[1:(n/2 + 1)])^2
}
