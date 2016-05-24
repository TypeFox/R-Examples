require(MASS)
require(lattice)
attach(mammals)
## trellis.device("pdf", color=FALSE, file = "bbMargins.pdf", width = 4,
##                height = 4)
source("latticeMono.R") ## until it's in the package
invisible(latticeMono(xyplot(log(brain) ~ log(body)), gray(.5)))
## dev.off()

