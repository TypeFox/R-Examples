\dontrun{#REX
library(psd)

##
## Spline gradient
##

set.seed(1234)
x <- seq(0,5*pi,by=pi/64)
y <- cos(x) #**2

splineGrad(x, y, TRUE)

# unfortunately, the presence of
# noise will affect numerical derivatives
y <- y + rnorm(length(y), sd=.1)
splineGrad(x, y, TRUE)

# so change the smoothing used in smooth.spline
splineGrad(x, y, TRUE, spar=0.2)
splineGrad(x, y, TRUE, spar=0.6)
splineGrad(x, y, TRUE, spar=1.0)

}#REX
