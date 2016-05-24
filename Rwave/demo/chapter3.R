## Example 3.1 Gabor Functions
## Generate a Gabor function at given frequency (here 0.1 Hz), for a
## sampling frequency (here 1 Hz), location (here 256) and scale
## (here 25):
gab <- gabor(512, 256, .1, 25)
par(mfrow=c(1,2))
plot(Re(gab), type="l", ylab="")
title("Real part")
plot(Im(gab), type="l", ylab="")
title("Imaginary part")

## Example 3.2 Transients
## Compute Gabor transform of transient signal, here with n_f = 50,
## delta_f = 0.02 and two different values for the scale parameter,
## namely scale = 10, 18:
par(mfrow=c(4,1))
data(A0)
plot(A0, type="l", xaxs="i")
title("Transients")
cgtA0 <- cgt(A0, 50, .02, 10)
cgtA0 <- cgt(A0, 50, .02, 18)
## To display the phase of the Gabor transform
tmp <- cleanph(cgtA0, .1)
title("Phase of the Gabor transform")

## Example 3.3 Sine wave
## Generate the sequence, here a sine wave of frequency 1/32 Hz sampled
## with unit sampling frequency
x <- 1:512
sinewave <- sin(2*pi*x/32)
par(mfrow=c(3,1))
plot(sinewave, type="l")
title("Sine wave")
## Compute the Gabor transform with n_f = 50, delta_f = .005 and scale
## sigma = 25.  This corresponds to frequencies ranging from 0 to 0.125
## Hz.  Display the phase:
cgtsinewave <- cgt(sinewave, 50, .005, 25)
tmp <- cleanph(cgtsinewave, .01)
title("Gabor Transform Phase")

## Example 3.4 Chirp
## Generate the chirp and compute the Gabor transform between frequencies
## 0 and 0.125 Hz:
x <- 1:512
chirp <- sin(2*pi*(x + 0.002*(x-256)*(x-256))/16) 
par(mfrow=c(3,1))
plot(ts(chirp), xaxs="i", xlab="", ylab="")
title('Chirp signal')
## The result is displayed in Figure 3.5

mycgt <- function (input, nvoice, freqstep = (1/nvoice),
                   scale = 1, plot = TRUE) {
  oldinput <- input
  isize <- length(oldinput)
  tmp <- adjust.length(oldinput)
  input <- tmp$signal
  newsize <- length(input)
  pp <- nvoice
  Routput <- matrix(0, newsize, pp)
  Ioutput <- matrix(0, newsize, pp)
  output <- matrix(0, newsize, pp)
  dim(Routput) <- c(pp * newsize, 1)
  dim(Ioutput) <- c(pp * newsize, 1)
  dim(input) <- c(newsize, 1)
  z <- .C("Sgabor", as.double(input), Rtmp = as.double(Routput), 
          Itmp = as.double(Ioutput), as.integer(nvoice), as.double(freqstep), 
          as.integer(newsize), as.double(scale), PACKAGE="Rwave")
  Routput <- z$Rtmp
  Ioutput <- z$Itmp
  dim(Routput) <- c(newsize, pp)
  dim(Ioutput) <- c(newsize, pp)
  i <- sqrt(as.complex(-1))
  output <- Routput[1:isize, ] + Ioutput[1:isize, ] * i
  if (plot) {
    image(1:newsize, seq(0, nvoice*freqstep/2, length=nvoice),
          Mod(output), xlab = "Time", ylab = "Frequency")
    title("Gabor Transform Modulus")
  }
  output
}

cgtchirp <- mycgt(chirp, 50, .005, 25)
tmp <- cleanph(cgtchirp, .01)
