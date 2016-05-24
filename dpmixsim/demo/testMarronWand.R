## Testing dpmixsim discriminating power for
## Marron Wand Densities as defined in Maechler's nor1mix

require(nor1mix)
op <- par(no.readonly = TRUE, ask=FALSE)
## These are defined as norMix() calls in  ../R/zMarrWand-dens.R
nms <- ls(pat="^MW.nm", "package:nor1mix")
options(warn=-1)
nms <- nms[order(as.numeric(substring(nms,6)))]
options(warn=0)
## Plot all of them:
par(mfrow=c(4,4), mgp = c(1.2, 0.5, 0), tcl = -0.2,
          mar = .1 + c(2,2,2,1), oma = c(0,0,3,0))
for(n in nms[-17]) plot(get(n, "package:nor1mix"))
mtext("The Marron-Wand Densities", outer= TRUE, font= 2, cex= 1.6)
par(op)

## test 
require(dpmixsim)
n <- as.numeric(readline("Choose one MW density to simulate {1:16}, or other to quit: "))
kmax   <- 100
while(n %in% 1:16) {
  mw <- get(nms[n], "package:nor1mix")
  if(n < 3) {
    b <- 0.1
    minvar <- 0.005
  }
  else {
	  if(n < 10) {
 	   b <- 0.01
 	   minvar <- 0.001
 	  }
    else {
 	    b <- 0.001
 	    minvar <- 0.0005
 	  }
	}
  x <- rnorMix(2000, mw)
  x <- prescale(x) 
  cat("simulating density MW",n,"\n")
  res <- dpmixsim(x, M=1, a=1, b=b, upalpha=1, maxiter=3000, rec=2000, kmax=kmax, nclinit=10, minvar=minvar)
  if(length(dev.list()) == 1) {dev.new(); par(ask=FALSE)}
  z <-  postdpmixciz(x, res=res, rec=2000, ngrid=200, plot=T, kmax=kmax)
  ##
  n <- as.numeric(readline("Choose one MW density to test {1:16}, or other to quit: "))
}
graphics.off()

