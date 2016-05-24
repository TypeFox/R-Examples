## Copyright (C) 1999-2000 Paul Kienzle
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## usage: [S [, f [, t]]] = specgram(x [, n [, Fs [, window [, overlap]]]])
##
## Generate a spectrogram for the signal. This chops the signal into
## overlapping slices, windows each slice and applies a Fourier
## transform to determine the frequency components at that slice.
##
## x: vector of samples
## n: size of fourier transform window, or [] for default=256
## Fs: sample rate, or [] for default=2 Hz
## window: shape of the fourier transform window, or [] for default=hanning(n)
##    Note: window length can be specified instead, in which case
##    window=hanning(length)
## overlap: overlap with previous window, or [] for default=length(window)/2
##
## Return values
##    S is complex output of the FFT, one row per slice
##    f is the frequency indices corresponding to the rows of S.
##    t is the time indices corresponding to the columns of S.
##    If no return value is requested, the spectrogram is displayed instead.
##
## Example
##    x = chirp([0:0.001:2],0,2,500);  # freq. sweep from 0-500 over 2 sec.
##    Fs=1000;                  # sampled every 0.001 sec so rate is 1 kHz
##    step=ceil(20*Fs/1000);    # one spectral slice every 20 ms
##    window=ceil(100*Fs/1000); # 100 ms data window
##    specgram(x, 2^nextpow2(window), Fs, window, window-step)
##
##    ## Speech spectrogram
##    [x, Fs] = auload(file_in_loadpath("sample.wav")); # audio file
##    step = fix(5*Fs/1000);     # one spectral slice every 5 ms
##    window = fix(40*Fs/1000);  # 40 ms data window
##    fftn = 2^nextpow2(window); # next highest power of 2
##    [S, f, t] = specgram(x, fftn, Fs, window, window-step)
##    S = abs(S(2:fftn*4000/Fs,:)); # magnitude in range 0<f<=4000 Hz.
##    S = S/max(S(:));           # normalize magnitude so that max is 0 dB.
##    S = max(S, 10^(-40/10));   # clip below -40 dB.
##    S = min(S, 10^(-3/10));    # clip above -3 dB.
##    imagesc(flipud(log(S)));   # display in log scale
##
##
##     imagesc(flipud(log(S(idx,:))))

## 2001-07-05 Paul Kienzle <pkienzle@users.sf.net>
## * remove "See also spectrogram"
## * add notes on selecting parameters for the spectrogram

specgram <- function(x, n = min(256, length(x)), Fs = 2, window = hanning(n),
     overlap = ceiling(length(window)/2)){ 
  ## if only the window length is given, generate hanning window

  if(!is.numeric(x))
    stop("'x' has to be a numeric.")
      
  if (length(window) == 1)
    window <- hanning(window)

  ## should be extended to accept a vector of frequencies at which to
  ## evaluate the fourier transform (via filterbank or chirp
  ## z-transform)
  if (length(n) > 1)
    stop("specgram does not handle frequency vectors yet")

  ## compute window offsets
  win_size <- length(window)
  if (win_size > n) {
    n <- win_size
    warning("specgram fft size adjusted to", n)
  } 
  step <- win_size - overlap

  ## build matrix of windowed data slices
  if (length(x) > win_size)
    offset <- seq(1, length(x)-win_size, by = step)
  else
    offset <- 1    
  S <- matrix(0, n, length(offset))
  for (i in seq_along(offset)) {
    S[1:win_size, i] <- x[offset[i]:(offset[i] + win_size - 1)] * window
  }

  ## compute fourier transform
  S <- mvfft(S)

  ## extract the positive frequency components
  if (n %% 2 == 1)
    ret_n <- (n+1)/2
  else
    ret_n <- n/2
  S <- S[1:ret_n, ]
  
  f <- (0:(ret_n-1))*Fs/n
  t <- offset/Fs
  
  res <- list(S = S, f = f, t = t)
  class(res) <- "specgram"
  res
}

print.specgram <- plot.specgram <- function(x, col = gray(0:512 / 512), xlab="time", ylab="frequency", ...) {
  image(x$t, x$f, 20 * log10(t(abs(x$S))), col = col, xlab = xlab, ylab = ylab, ...)
}


###specgram(chirp(seq(-2, 15, by = .001), 400, 10, 100, 'quadratic'))
###specgram(chirp(seq(0, 5, by = 1/8000), 200, 2, 500, "logarithmic"), Fs = 8000)
###
###
###Fs=1000
###x = chirp(seq(0,2, by = 1/Fs), 0, 2, 500)   # freq. sweep from 0-500 over 2 sec.
###step = ceiling(20*Fs/1000)      # one spectral slice every 20 ms
###window = ceiling(100*Fs/1000)   # 100 ms data window
###specgram(x)

#! ## test of returned shape
#!assert (rows(S), 128)
#!assert (columns(f), rows(S))
#!assert (columns(t), columns(S))
#!test [S, f, t] = specgram(x')
#!assert (rows(S), 128)
#!assert (columns(f), rows(S))
#!assert (columns(t), columns(S))
#!error (isempty(specgram([])))
#!error (isempty(specgram([1, 2 ; 3, 4])))
#!error (specgram)

#!demo
#! Fs=1000
#! x = chirp([0:1/Fs:2],0,2,500);  # freq. sweep from 0-500 over 2 sec.
#! step=ceil(20*Fs/1000);    # one spectral slice every 20 ms
#! window=ceil(100*Fs/1000); # 100 ms data window
#!
#! ## test of automatic plot
#! [S, f, t] = specgram(x)
#! specgram(x, 2^nextpow2(window), Fs, window, window-step)
#! disp("shows a diagonal from bottom left to top right")
#! input("press enter:","s")
#!
#! ## test of returned values
#! S = specgram(x, 2^nextpow2(window), Fs, window, window-step)
#! imagesc(20*log10(flipud(abs(S))))
#! disp("same again, but this time using returned value")

#!demo
#! ## Speech spectrogram
#! [x, Fs] = auload(file_in_loadpath("sample.wav")); # audio file
###data(sampleWav)
###Fs = wav$rate
###step = trunc(5*Fs/1000);     # one spectral slice every 5 ms
###window = trunc(40*Fs/1000);  # 40 ms data window
###fftn = 2^nextpow2(window); # next highest power of 2
###spg = specgram(wav$sound, fftn, Fs, window, window-step)
###S = abs(spg$S[2:(fftn*4000/Fs),]) # magnitude in range 0<f<=4000 Hz.
###S = S/max(S);         # normalize magnitude so that max is 0 dB.
###S[S < 10^(-40/10)] = 10^(-40/10)   # clip below -40 dB.
###S[S > 10^(-3/10)] = 10^(-3/10))    # clip above -3 dB.
###image(t(20*log10(S)), axes = FALSE) #, col = gray(0:255 / 255))
#!
#! % The image contains a spectrogram of 'sample.wav'
