###
### $Id: denoise.R 58 2014-06-24 20:27:34Z plroebuck $
### Package Demo
###

library("rwt")

message("************************************")
message("Creating doppler signal for test")

sig <- makesig(SIGNAL.DOPPLER)
s <- as.vector(sig$x)
N <- sig$N
n <- rnorm(N)

x <- s + n/10     # (approximately 10dB SNR)

message("Get Daubechies's scaling filter")
h <- daubcqf(8)$h.0

## Denoise 'x' with the default method based on the DWT
message("Denoise signal using DWT")
#debug(.dwt)
#debug(mdwt)
#debug(midwt)
ret.dwt <- denoise.dwt(x, h)
xd <- ret.dwt$xd
xn <- ret.dwt$xn
opt1 <- ret.dwt$option

## Denoise 'x' using the undecimated (LSI) wavelet transform
message("Denoise signal using UDWT")
#debug(.udwt)
#debug(mrdwt)
#debug(mirdwt)
ret.udwt <- denoise.udwt(x, h)
yd <- ret.udwt$xd
yn <- ret.udwt$xn
opt2 <- ret.udwt$option

## Plot results
message("Plot results")
dev.new()
plotSignalTransformation(x, s, "Original")
dev.new()
plotSignalTransformation(xd, s, "Decimated Wavelet Transform")
dev.new()
plotSignalTransformation(yd, s, "Undecimated Wavelet Transform")

