## ----eval=TRUE, echo=TRUE, label="Filtered signal analogy."--------------
set.seed(1234)
N <- 1028
x <- rnorm(N, mean = 0, sd = 1)
# Load signal processing
library(signal, warn.conflicts=FALSE)
# Construct an FIR filter
f <- c(0, 0.2, 0.2, 0.3, 0.3, 0.5)*2
m <- c(0, 0, 1, 1, 0, 0)
fw <- signal::fir2(N, f, m)
# complex filter response
fh <- signal::freqz(fw, Fs=1)
f <- c(0, 0.12, 0.12, 0.22, 0.22, 0.5)*2
fwl <- signal::fir2(N, f, m)
fhl <- signal::freqz(fwl, Fs=1)
f <- c(0, 0.28, 0.28, 0.38, 0.38, 0.5)*2
fwu <- signal::fir2(N, f, m)
fhu <- signal::freqz(fwu, Fs=1)
# convolution
xf <- signal::filter(fw, x)
# PSD using stats::spectrum
Sx <- spectrum(x, pad=1, plot=FALSE, taper=0.2)
Sxf <- spectrum(xf, pad=1, plot=FALSE, taper=0.2)

## ----eval=TRUE, echo=FALSE, fig.width=6, fig.height=3, label=FILTERS-----
Sx$df <- Sx$bandwidth <- NA
par(mar=c(0,0,2,0), oma=rep(0,4))
plot(Sx, col="grey", log="dB", ylim=c(-150,20), xlim=c(0.1,0.4),
xlab="", ylab="", ci.col=NA, axes=FALSE,
main="PSDs of raw and filtered process")
lines(fhu$f, 20*log10(Mod(fhu$h)), col="red", lty=3, lwd=2)
lines(fhl$f, 20*log10(Mod(fhl$h)), col="red", lty=3, lwd=2)
lines(fh$f, 20*log10(Mod(fh$h)), col="red", lwd=2)
plot(Sxf, log="dB", add=TRUE)

## ----eval=TRUE, echo=TRUE, label="Synthetic white noise and a DFT."------
# using x from the filter analogy section
xv <- var(x)
X <- fft(x)
class(X)
length(X)

## ----eval=TRUE, echo=TRUE, label="Amplitude and phase spectra."----------
Sa <- Mod(X) # Amplitude spectrum
Sp <- Arg(X) # Phase spectrum
XC <- Conj(X)
all.equal(Se <- Sa**2, Se_2 <- Mod(XC * X), Se_2R <- Mod(X * XC))

## ----eval=TRUE, echo=TRUE, label="Nyquist frequencies."------------------
fsamp <- 1  # sampling freq, e.g. Hz
fNyq <- fsamp/2   # Nyquist freq
Nf <- N/2  # number of freqs
nyfreqs <- seq.int(from=0, to=fNyq, length.out=Nf)
S <- Se[2:(Nf+1)] * 2 / N   # Finally, the PSD!

## ----eval=TRUE, echo=TRUE, label=OPTIM_EXPECT----------------------------
# 0) Setup optimization function for dof, using conjugate gradients\\
#    min L1 |PSD - Chi^2(dof)|
Chifit <- function(PSD){optim(list(df=mean(PSD)), function(dof){
  sum(log(PSD)) - sum(log(dchisq(PSD, dof))) }, method="CG")}
# 1) run optimization
Schi <- Chifit(S)
# Get 'df', the degrees of freedom
print(dof <- Schi$par[[1]]) 

## ----eval=TRUE, echo=TRUE, label=MEAN_EXPECT-----------------------------
# compare with the mean and median
c(mSn <- mean(S), median(S))

## ----eval=TRUE, echo=FALSE, fig.width=5, fig.height=5, label=QQFIT-------
par(pty="s", las=1)
ttl <- expression("Q-Q plot for PSD and" ~~ chi^2 ~~ "fit")
qqplot(log(qchisq(ppoints(N), df=dof)), log(S), main = ttl, ylab="Spectrum quantiles", xlab="Distribution quantiles")
abline(0,1, col="red")

## ----eval=TRUE, echo=FALSE, fig.width=6, fig.height=3.5, label=PSD-------
par(las=1, mgp = c(2.2, 1, 0))
ylab <- expression("units"^2 ~~ "/ frequency")
plot(nyfreqs, S, type="h", xlab="Nyquist frequency", ylab=ylab, yaxs="i")
abline(h=dof, lwd=2, col="red")

## ----eval=TRUE, echo=TRUE, label="Test normalization."-------------------
mSn <- dof
test_norm <- function(sval, nyq, xvar){svar <- sval * nyq; return(svar/xvar)}
print(xv_1 <- test_norm(mSn, fNyq, xv))
xv_2 <- sum(S)/Nf * fNyq / xv  # an alternate test
all.equal(xv_1, xv_2)

## ----eval=TRUE, echo=TRUE, label="Apply correct normalization."----------
fsamp <- 20
fNyq <- fsamp / 2
freqs <- fsamp * nyfreqs 
Snew <- S / fsamp
# Test variance crudely
mSn <- mean(Snew)
test_norm(mSn, fNyq, xv)

## ----eval=TRUE, echo=TRUE, label="DB"------------------------------------
# decibel function
dB <- function(y) 10*log10(y)

## ----eval=TRUE, echo=FALSE, fig.width=6, fig.height=3.8, label=PSD2------
par(las=1)
plot(freqs, dB(S), type="l", col="dark grey", xlab="Frequency", ylab="dB")
lines(freqs, dB(Snew), lwd=1.3)
lines(c(0,fNyq), rep(dB(mSn),2), lwd=3, col="red")
abline(h=dB(1/fNyq), col="blue")

## ----eval=TRUE, echo=TRUE, label="Change the sampling frequency."--------
fsamp <- 20
xt <- ts(x, frequency=fsamp)
pgram20 <- spec.pgram(xt, pad=1, taper=0, plot=FALSE)
pgram01 <- spec.pgram(ts(xt, frequency=1), pad=1, taper=0, plot=FALSE)

## ----eval=TRUE, echo=TRUE------------------------------------------------
mSn/mean(pgram20$spec)

## ----eval=TRUE, echo=FALSE, fig.width=6, fig.height=3.5, label=NORMS-----
par(las=1)
plot(pgram01, log="dB", xlim=c(0,10), ylim=36*c(-1,.3), main="", col="dark grey")
plot(pgram20, log="dB", add=TRUE)
abline(h=-dB(c(1, 20)*2), col=c("dark grey","black"))
abline(v=.5*c(1,20), lty=3)
#lines(c(0,fNyq), rep(dB(mSn),2), lwd=1.5, col="red")
abline(h=dB(1/fNyq), col="blue")

## ----eval=TRUE, echo=TRUE, label="Test the normalization again."---------
test_norm(mean(pgram01$spec), 0.5, xv)
test_norm(mean(pgram20$spec), 10, xv)

## ----eval=TRUE, echo=TRUE, label="Double sided PSD from spectrum."-------
psd1 <- spec.pgram(x, plot=FALSE)
psd2 <- spec.pgram(xc<-complex(real=x, imag=x), plot=FALSE, demean=TRUE)
mx <- mean(Mod(x))
mxc <- mean(Mod(xc))
(mxc/mx)**2
mean(psd2$spec / psd1$spec)

## ----eval=TRUE, echo=TRUE, label=SI--------------------------------------
utils::sessionInfo()

