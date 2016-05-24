#Deconvolution of  monthly sunspots data with Lorentzian profile
library(datasets)
library(Peaks)
library(stats)
library(graphics)
#Generate spectrum
abs(fft(spec.taper(as.vector(sunspot.month),p=0.5)))->smf
#Remove background
p <- smf-SpectrumBackground(smf,iterations=100)
#Generate Lorentzian responce vector with FWHM=pi*5
x <- 1:100
Z<-1/(1+((x-pi*5)/5)^2)
#Deconvolve low-frequency region of spectrum
l <- SpectrumDeconvolution(p[1:500],Z,iterations=20,boost=1.1,repetitions=4)
plot(p,type="p",xlim=c(0,200))
lines(l,col="red")
