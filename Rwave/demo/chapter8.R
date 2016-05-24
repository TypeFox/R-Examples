## Example 8.1 "Pure" Transient Analysis/Synthesis
## Compute the dyadic wavelet transform and the corresponding local
## extrema:
data("C0")
dwC0 <- mw(C0,5)
extC0 <- ext(dwC0)
## Reconstruction from the extrema:
recC0 <- mrecons(extC0)
## Example 8.2 "Noisy" Transient Analysis
data("C4")
dwC4 <- mw(C4,5)
extC4 <- ext(dwC4)
## Example 8.3 : "Noisy" Transient Detection/Synthesis
## Trim the extrema:
trC4 <- mntrim(extC4)
## Reconstruction from the trimmed extrema:
recC4 <- mrecons(trC4)
## Example 8.5 Simple Reconstruction of a Speech Signal
data("HOWAREYOU")
plot.ts(HOWAREYOU)
cgtHOW <- cgt(HOWAREYOU,70,0.01,60)
clHOW <- crc(Mod(cgtHOW),nbclimb=2000)
## Simple reconstruction:
srecHOW <- scrcrec(HOWAREYOU,cgtHOW,clHOW,ptile=0.0001,plot=2)
