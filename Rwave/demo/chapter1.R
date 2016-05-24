## Example 1.3 Dolphin click
## Read and plot the dolphin click data
data("click")
fclick <- fft(click)
par(mfrow=c(2,1))
plot.ts(click)
plot.ts(Mod(fclick))
## Subsampling the data:
k <- 1:1249
sclick <- click[2*k]
plot.ts(click)
plot.ts(sclick)
## Fourier transforms of the original and click data:
fsclick <- fft(sclick)
nfclick <- c(fclick[1250:2499], fclick[1:1249])
nfsclick <- c(fsclick[625:1249], fsclick[1:624])
plot.ts(Mod(nfclick), xaxt="n")
axis(1, at=c(0,500,1000,1500,2000,2500),
     labels=c("-1250","-750","-250","250","750","1250"))
plot.ts(Mod(nfsclick), xaxt="n")
axis(1, at=c(0,200,400,600,800,1000,1200),
     labels=c("-600","-400","-200","0","200","400","600"))
par(mfrow=c(1,1))
