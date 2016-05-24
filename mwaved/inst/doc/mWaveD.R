## ----setup, echo = FALSE, cache = TRUE-----------------------------------
library(ggplot2)

shift <- function(x) {
  x <- as.matrix(x)
  n2 <- dim(x)[1]/2 + 1
  y <- rbind(x[(n2+1):n,], x[1:n2,])
  y
}

shiftdom <- function(x) {
  n2 <- length(x)/2 + 1
  y <- c(x[(n2+1):n], x[1:n2])
  y
}

l <- labs(x = 'x', y = '')
blank <-  theme(legend.position = "none")

ggplotTime <- function(G, x = 'time', title, alpha = 1) {
  dat <- switch(x, 
                time = stack(as.data.frame(G)),
                shiftedTime  = stack(as.data.frame(shift(G))),
                freq = stack(as.data.frame(Mod(G))),
                shiftedFreq = stack(as.data.frame(shift(Mod(G)))))
  m <- ncol(G)
  n <- nrow(G)
  dat$x = rep(switch(x, 
                     time = seq(from = 0, to = 1 - 1/n, length = n),
                     shiftedTime = seq(from = -0.5 + 1/n, to = 0.5, length = n),
                     freq = 0:(n-1),
                     shiftedFreq = -(floor(n/2) - 1):floor(n/2)), m)
  p <- ggplot(dat) + geom_line(aes(x = x, y = values, colour = ind), alpha = alpha) + ggtitle(title) + l + blank
  if (x == 'shiftedFreq') {
    p <- p + xlim(-100,100)
  }
    
  p
}


## ----firsteg, fig.show = 'hide'------------------------------------------
library(mwaved)
n <- 1024
m <- 3
signal <- makeDoppler(n)
x <- seq(from = 0, to = 1 - 1/n, length = n)
Gbox <- boxcarBlur(n, width = 1/sqrt(c(13, 343, 89)))
Xbox <- blurSignal(signal, Gbox)
head(Xbox)
head(Gbox)
plot(x, signal, type = 'l', main  = 'Doppler signal, f')
matplot(x, Xbox, type = 'l', main = 'g[l] * f')
matplot(x, Gbox, type = 'l', main = 'g[l]')
shift <- function(x) {
  x <- as.matrix(x)
  n2 <- dim(x)[1]/2 + 1
  y <- rbind(x[(n2+1):n,], x[1:n2,])
  y
}
x <- seq(from = -0.5 - 1/n, to = 0.5, length = n)
matplot(x, shift(Gbox), type = 'l', main = 'Shifted g[j]')

## ----dataplots, fig.show = 'hold', echo = FALSE, results = 'hide'--------
library(mwaved)
library(ggplot2)
n <- 1024
m <- 3
shift <- function(x) {
  x <- as.matrix(x)
  n2 <- dim(x)[1]/2 + 1
  y <- rbind(x[(n2+1):n,], x[1:n2,])
  y
}
signal <- makeDoppler(n)
l <- labs(x = 'x', y = "") 
x <- seq(from = 0, to = 1 - 1/n, length = n)
Gbox <- boxcarBlur(n, width = 1/sqrt(c(13, 343, 89)))
Xbox <- blurSignal(signal, Gbox)
ggplot(data.frame(x = x, y = signal)) + geom_line(aes(x = x, y = y)) + ggtitle("Doppler signal, f") + l
ggplotTime(Xbox, title = "g[l] * f")
ggplotTime(Gbox, title = "g[l]")
ggplotTime(Gbox, x = 'shiftedTime', title = 'Shifted g[l]')

## ----simulation-functions, fig.show='hide'-------------------------------
library(mwaved)
n <- 1024
x <- seq(from = 0, to = 1 - 1/n, length = n)
Doppler <- makeDoppler(n)
plot(x, Doppler, type = 'l')

## ----plot-simulation, echo = FALSE, results='hide', message = FALSE, fig.show='hold'----
library(mwaved)
library(ggplot2)
n <- 1024
dat <- data.frame(x = seq(from = 0, to = 1-1/n, length = n), 
                  Doppler = makeDoppler(n),
                  LIDAR = makeLIDAR(n),
                  Bumps = makeBumps(n),
                  Blocks = makeBlocks(n))
signals <- c("Doppler", "LIDAR", "Bumps", "Blocks")
baseP <- ggplot(dat)
plots <- lapply(1:4, function(i) baseP + geom_line(aes_string(x = 'x', y = signals[i])) + ggtitle(signals[i]) + ylab(""))
# do.call(grid.arrange,  plots)
plots

## ----blur-functions, fig.show='hide'-------------------------------------
library(mwaved)
n <- 1024
m <- 3
Gbox <- boxcarBlur(n, width = 1/sqrt(c(13, 343, 89)))
Gdir <- directBlur(n, m)
Gsmooth <- gammaBlur(n, shape = c(0.25, 0.5, 0.75), scale = rep(0.25, m))
x <- seq(from = 0, to = 1 - 1/n, length = n)

n2 <- floor(n/2)
w = -(n2 - 1):n2
Gdirfft <- mvfft(Gdir)
Gboxfft <- mvfft(Gbox)
Gsmoothfft <- mvfft(Gsmooth)
shift <- function(x) {
  x <- as.matrix(x)
  n2 <- dim(x)[1]/2 + 1
  y <- rbind(x[(n2+1):n,], x[1:n2,])
  y
}
xlims = c(-100,100)
matplot(x, Gdir, type = 'l', main = 'Direct blur matrices')
matplot(w, shift(Mod(Gdirfft)), type = 'l', main = 'Fourier shifted', xlim = xlims)
matplot(x, Gsmooth, type = 'l', main = 'Smooth blur matrices')
matplot(w, shift(Mod(Gsmoothfft)), type = 'l', main = 'Fourier shifted', xlim = xlims)
matplot(x, Gbox, type = 'l', main = 'Boxcar blur matrices')
matplot(w, shift(Mod(Gboxfft)), type = 'l', main = 'Fourier shifted', xlim = xlims)


## ----blur-plots, fig.show='hold', warning = FALSE, echo = FALSE, results = 'hide'----
ggplotTime(Gdir, title = 'Direct blur matrices')
ggplotTime(Gdirfft, x = 'shiftedFreq', title = 'Fourier shifted')
ggplotTime(Gsmooth, title = 'Smooth blur matrices')
ggplotTime(Gsmoothfft, x = 'shiftedFreq', title = 'Fourier shifted')
ggplotTime(Gbox, title = 'Boxcar blur matrices')
ggplotTime(Gboxfft, x = 'shiftedFreq', title = 'Fourier shifted')


## ----noise-sim, fig.show='hide'------------------------------------------
library(mwaved)
n <- 1024
m <- 3
alpha <- seq(from = 0.5, to = 1, length = m)
x <- seq(from = 0, to = 1 - 1/n, length = n)
sigma <- c(1, 0.75, 0.5)
E <- multiNoise(n, sigma = sigma, alpha = alpha)
matplot(x, E, type = 'l', main = 'Long Memory Noise')
DE <- multiNoise(n, sigma = sigmaSNR(signal = Xbox, SNR = c(10, 20 ,30)), alpha = alpha)
Y <- Xbox + DE
matplot(x, Y, type = 'l', main = 'Y{i,l] = g[j]*f + e[i,l]')

## ----noise-gg-plots, echo = FALSE, message = FALSE, fig.show='hold'------
ggplotTime(E, title = 'Long Memory Noise', alpha = 0.5)
ggplotTime(Y, title = 'Y[i,l] = g[l]*f + e[i,l]')

## ----estimation-example, fig.show='hide'---------------------------------
library(mwaved)
n <- 1024
m <- 3
alpha <- c(0.95, 0.8, 0.75)
SNR <- c(30, 28, 25)
signal <- makeHeaviSine(n)
Gsmooth <- gammaBlur(n, shape = c(0.5, 0.66, 0.75), scale = rep(0.25, 3))
# G <- boxcarBlur(n)
X <- blurSignal(signal, Gsmooth)
E <- multiNoise(n, sigma = sigmaSNR(X, SNR))
Ysmooth <- X + E
x <- seq(from = 0, to = 1 - 1/n, length = n)
signalEstimate <- multiEstimate(Ysmooth, Gsmooth, alpha)
head(signalEstimate)
matplot(x, Ysmooth, type = 'l')
plot(x, signal, type = 'l')
lines(x, signalEstimate, col = 2)

## ----ggplot-multiestimate, fig.show = 'hold', echo = FALSE, results = 'hide'----
ggplotTime(Ysmooth, title = 'Multichannel signal')
ggplotTime(cbind(signal, signalEstimate), title = "TI mWaveD estimate")

## ----full-multiWaveD, fig.show = 'hide'----------------------------------
signalmWaveD <- multiWaveD(Ysmooth, Gsmooth, alpha)
names(signalmWaveD)
plot(signalmWaveD)

## ----ggplot-multiWaveD, fig.show = 'hold', results = 'hide', echo = FALSE----
plot(signalmWaveD, which = 1)
plot(signalmWaveD, which = 2)
plot(signalmWaveD, which = 3)
plot(signalmWaveD, which = 4)

## ----single-multiWaveD-plot----------------------------------------------
plot(signalmWaveD, which = 4)

## ----summary-mWaveD, results = 'markup'----------------------------------
summary(signalmWaveD)

## ----smooth-plots, fig.show = 'hold'-------------------------------------
smoothObj <- multiWaveD(Ysmooth, Gsmooth, alpha, resolution = 'smooth')
plot(smoothObj, which = 2)
plot(smoothObj, which = 3)

Xbox <- blurSignal(signal, Gbox)
E <- multiNoise(n, sigma = sigmaSNR(Xbox, SNR))
Ybox <- Xbox + E

boxObjSmooth <- multiWaveD(Ybox, Gbox, alpha = 1, resolution = 'smooth')
bsj1 <- boxObjSmooth$j1
boxObjBlock <- multiWaveD(Ybox, Gbox, alpha = 1, resolution = 'block')
bbj1 <- boxObjBlock$j1
paste0('estimate j1 = ',bsj1, ' or ',bbj1, ' for the smooth and block method respectively.')
plot(boxObjSmooth, which = 2)
plot(boxObjSmooth, which = 3)
plot(boxObjBlock, which = 2)
plot(boxObjBlock, which = 3)

## ----noisySmooth, fig.show = 'hold'--------------------------------------
NGsmooth <- Gsmooth + multiNoise(n, sigma = sigmaSNR(Gsmooth, rep(10,3)))
noisySmoothObj <- multiWaveD(Ysmooth, NGsmooth, alpha, resolution = 'smooth')
plot(smoothObj, which = 3)
plot(noisySmoothObj, which = 3)

## ----blockPlots, fig.show = 'hold'---------------------------------------
plot(boxObjBlock, which = 2)
plot(boxObjBlock, which = 3)

## ----thresholding-plots, results='hide', echo = FALSE--------------------
  library(ggplot2)
  t = 3
  thr = 1
  l <- labs(x = 'x', y = '') 
  xt <- seq(-t, t, length = 2^10)
  yt <- rep(0, length(xt))
  ht <- xt*(abs(xt) > thr)
  st <- (abs(xt) > thr)*(xt - (xt > thr) + (xt < -thr))
  gt <- (xt - thr^2/xt)*(abs(xt) > thr)
  threshDat <- data.frame(x = xt, ht = ht, st = st, gt = gt)
  baseT <- ggplot(threshDat) + l
  
  hardPlot <- baseT + geom_line(aes(x = x, y = ht)) + geom_line(aes(x = x, y = x), lty = 'dashed') + ggtitle('Hard thresholding')
  softPlot <- baseT + geom_line(aes(x = x, y = st)) + geom_line(aes(x = x, y = x), lty = 'dashed') + ggtitle('Soft thresholding')
  garrotePlot <- baseT + geom_line(aes(x = x, y = gt)) + geom_line(aes(x = x, y = x), lty = 'dashed') + ggtitle('Garrote thresholding')
hardPlot
softPlot
garrotePlot

## ----waveletCoef-eg------------------------------------------------------
wave <- multiCoef(Ysmooth, Gsmooth, alpha = alpha)
str(wave)
names(wave)

## ----plot-waveCoef, fig.show = 'hold'------------------------------------
library(mwaved)
n <- 1024
signal <- makeLIDAR(n)
Y <- signal + multiNoise(n, sigma = sigmaSNR(signal, 20))
# Without thresholding
rawWave <- multiCoef(Y)
thresh <- multiThresh(Y)
# With thresholding
threshWave <- waveletThresh(rawWave, thresh)
sigCoefs <- multiCoef(signal)
plot(sigCoefs)
plot(rawWave, threshWave)

## ----multiProj, fig.show='hide'------------------------------------------
x <- seq(from = 0, to = 1 - 1/n, length = n)
plot(x, multiProj(sigCoefs), ty = 'l', main = 'signal coefs')
plot(x, multiProj(rawWave), ty = 'l', main = 'noisy signal coefs')
plot(x, multiProj(threshWave), ty = 'l', main = 'thresholded signal coefs')
plot(x, multiEstimate(Y), ty = 'l', main = 'TI mWaveD estimate')

## ----multiProj-gg, fig.show='hold', echo = FALSE, results='hide'---------
ggplotTime(as.matrix(multiProj(sigCoefs)), title = 'signal coefs')
ggplotTime(as.matrix(multiProj(rawWave)), title = 'noisy signal coefs')
ggplotTime(as.matrix(multiProj(threshWave)), title = 'thresholded signal coefs')
ggplotTime(as.matrix(multiEstimate(Y)), title = 'TI mWaveD estimate')

## ----multiSigma----------------------------------------------------------
library(mwaved)
sigma <- c(0.5, 0.75, 2)
n <- 1024
E <- multiNoise(n, sigma = sigma)
multiSigma(E)

