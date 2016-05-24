require(seewave); data(tico); data(orni); data(pellucens)

op1<-par(ask=TRUE)

# different oscillograms of a tropical sparrow song
oscillo(tico,f=22050)
oscillo(tico,f=22050,k=2,j=2,byrow=TRUE)
op<-par(bg="grey")
oscillo(tico,f=22050,k=4,j=1,title=TRUE,colwave="black",
    coltitle="yellow",collab="red",colline="white",
    colaxis="blue",coly0="grey50")
par(op)

# overplot of oscillographic and envelope representation
oscillo(tico,f=22050)
par(new=TRUE)
env(tico,f=22050,colwave=2)

# temporal automatic measurements
timer(orni,f=22050,threshold=5,msmooth=c(40,0),
        bty="l",colval="blue")
title(main="Timer() for automatic time measurements",col="blue")

# instantaneous frequency
ifreq(tico,f=22050,threshold=5)
title(main="Instantaneous frequency using Hilbert transform")

# comparaison of a full spectrum and a mean spectrum of a cicada song
op<-par(mfrow=c(2,1))
spec(orni,f=22050,type="l")
title("spec()")
meanspec(orni,f=22050,wl=512,type="l")
title("meanspec()")
par(op)

# basic 2D spectrogram of a bird song
op <- par(op)
spectro(tico,f=22050,wl=512,ovlp=50,zp=16,collevels=seq(-40,0,0.5))
par(op)

# spectrogram and dominant frequency overlaid of a bird song
op <- par(op)
spectro(tico, f=22050, ovlp=50, palette=reverse.gray.colors.2, scale=FALSE)
par(new=T)
dfreq(tico, f=22050, ovlp=50, threshold=6, col="red", ann=FALSE, xaxs="i", yaxs="i")
par(op)

# 2D spectrogram of a cricket song with colour modifications
op <- par(op)
pellu2<-cutw(pellucens,f=22050,from=1,plot=FALSE)
spectro(pellu2,f=22050,wl=512,ovlp=85,collevels=seq(-25,0,1),osc=TRUE,palette=reverse.heat.colors,
colgrid="white", colwave="white",colaxis="white",collab="white", colbg="black")
par(op)

# sound synthesis
op <- par(op)
F1<-synth(f=22050,am=c(50,10),cf=2000,d=1,fm=c(500,5,0),plot=FALSE)
F2<-synth(f=22050,a=0.8,cf=4000,am=c(50,10),d=1,fm=c(500,5,0),plot=FALSE)
F3<-synth(f=22050,a=0.6,cf=6000,am=c(50,10),d=1,fm=c(500,5,2000),plot=FALSE)
F4<-synth(f=22050,a=0.4,cf=8000,am=c(50,10),d=1,fm=c(500,5,2000),plot=FALSE)
final1<-F1+F2+F3+F4
spectro(final1,f=22050,wl=512,ovlp=75,osc=TRUE)
title(main="synthesis of a AM/FM sound")
par(op)

# 3D spectrogram of a tropical sparrow  song
spectro3D(tico,f=22050,wl=512,ovlp=75,zp=16,maga=2)
