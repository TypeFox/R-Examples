\dontrun{#REX
library(psd)

##
## Riedel-Sidorenko-Parker taper optimization
##

set.seed(1234)
# some params
nd <- 512 # num data
ntap <- 10 # num tapers
nrm <- 40 # sharpness of the peaks rel 2*variance
#
# create a pseudo spectrum
# with broad peaks
x <- 0:(nd-1)
riex <- rnorm(nd) + nrm*abs(cos(pi*x/180) + 1.2)
riex <- riex + 8*nrm*dcauchy(x, nd/3)
riex <- riex + 5*nrm*dnorm(x, nd/2)
# some flat regions
riex[riex<25] <- 25
ried <- dB(riex, invert=TRUE)

# optimize tapers
rtap <- riedsid(riex, ntaper=ntap)
rtap2 <- riedsid2(riex, ntaper=ntap)

# plot
op <- par(no.readonly = TRUE)
par(mfrow=c(2,1), mar=rep(1.3,4), mai=rep(0.6,4))
# ... the mock spectrum
plot(riex, type="h", xaxs="i", ylim=c(0,200), main='Pseudo-spectrum') 
# ... the optimal tapers
plot(rtap2, main='Optimal tapers') 
# original tapers:
lines(rtap, col="red")
par(op)

}#REX
