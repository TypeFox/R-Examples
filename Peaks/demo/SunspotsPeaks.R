#Process monthly sunspots data
library(datasets)
library(Peaks)
library(stats)
library(graphics)
abs(fft(spec.taper(as.vector(sunspot.month),p=0.5)))->smf
SpectrumBackground(smf,iterations=100)->smb
plot(smf-smb,type="l",xlim=c(0,200))
SpectrumSearch(smf-smb)->z
lines(z$y,type="l",col="red")
points(y=rep(-10,length(z$pos)),x=z$pos,col="green",pch="+",cex=2)
print(paste(length(z$pos)/2," harmonics were found with ",z$pos[1],"-month base period",sep=""))
