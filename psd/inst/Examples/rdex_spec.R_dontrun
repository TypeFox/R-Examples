\dontrun{#REX
library(psd)

##
## Objects with class 'spec'
##

set.seed(1234)
xn <- rnorm(10)
x <- spectrum(xn, plot=FALSE)
xc <- psdcore(xn)

xdf <- as.data.frame(x)
str(xdf)
is.tapers(xdf$taper)

xdfc <- as.data.frame(xc)
str(xdfc)
is.tapers(xdfc$taper)

}#REX
