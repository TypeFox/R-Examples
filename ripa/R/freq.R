## $Header: /database/repository/rimage/R/freq.R,v 1.1.2.3 2004/03/17 06:35:18 tomo Exp $
##
## Copyright (c) 2003 Nikon Systems Inc.
## For complete license terms see file LICENSE


lowpass <- function(img, radius=40) {
  ## computes a low pass
  ## image is just for giving dimensions
  get.lp <- function(img, radius) {
    h <- dim(img)[1]
    w <- dim(img)[2]
    matrix(.C("getLowPass",
              as.double(img),
              as.integer(w),
              as.integer(h),
              as.integer(radius),
              spec = double(h*w),
              PACKAGE="ripa"
              )$spec, nrow=h, ncol=w)
  }
  if ( !is.na(dim(img)[3]) ){
	img = rgb2grey(I)
  }
  h <- dim(img)[1]
  w <- dim(img)[2]
  p = planFFT(h*w) 
  fft <- matrix(FFT(as.vector(t(img)), plan=p),nrow=h,ncol=w)
  t1 <- rbind(fft,fft)
  t2 <- cbind(t1,t1)
  x1 <- dim(fft)[1]/2 + 1
  x2 <- 1.5*dim(fft)[1]
  y1 <- dim(fft)[2]/2 + 1
  y2 <- 1.5*dim(fft)[2]
  fft <- t2[ x1:x2 , y1:y2 ]
  lp <- get.lp(img, radius)

  filtered <- fft*lp
  ifft <- matrix(IFFT(as.vector(t(filtered)),p),nrow=h,ncol=w)
  out1 = imagematrix(abs(ifft))
#img = rgb2grey(I)
}

highpass <- function(img, radius=40) {
  ## computes high pass
  get.hp <- function(img, radius) {
    h <- dim(img)[1]
    w <- dim(img)[2]
    matrix(.C("getHighPass",
              as.double(img),
              as.integer(w),
              as.integer(h),
              as.integer(radius),
              spec = double(w*h),
              PACKAGE="ripa"
              )$spec, nrow=dim(img)[1], ncol=dim(img)[2])
  }
  if (!is.na(dim(img)[3])){
	img = rgb2grey(I)
  }
  h <- dim(img)[1]
  w <- dim(img)[2]
  p = planFFT(h*w) 
  fft <- matrix(FFT(as.vector(t(img)), plan=p),nrow=h,ncol=w)
  t1 <- rbind(fft,fft)
  t2 <- cbind(t1,t1)
  x1 <- dim(fft)[1]/2 + 1
  x2 <- 1.5*dim(fft)[1]
  y1 <- dim(fft)[2]/2 + 1
  y2 <- 1.5*dim(fft)[2]
  fft <- t2[ x1:x2 , y1:y2 ]
  hp <- get.hp(img, radius)

  filtered <- fft*hp
  ifft <- matrix(IFFT(as.vector(t(filtered)),p),nrow=h,ncol=w)
  out1 = imagematrix(abs(ifft))
}


## the end of file
