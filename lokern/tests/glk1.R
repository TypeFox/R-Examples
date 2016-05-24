require(lokern)
data(xSim)
n <- length(xSim)
stopifnot(n == 75)
tt <- ((1:n) - 1/2)/n # equidistant x

(gk <- glkerns(tt, xSim))
summary(gk$est)
gk$bandwidth
glkerns(tt,xSim, deriv = 1)$bandwidth
glkerns(tt,xSim, deriv = 2)$bandwidth

demo("glk-derivs", ask = FALSE)# fails on Windows: grDevices::dev.interactive(TRUE)
##   ------------

p.3glks(tt, xSim, kord = 3)

p.3glks(tt, xSim, kord = 4, useB = 0.15)

str(p.3glks(tt, xSim, kord = 5, useB = 0.12) ) # k.ord = (4,5,4) => less sensiacl?

p.3glks(tt, xSim, kord = 6, useB = 0.2, derivs = 0:3) # k.ord = (6,5,6, 5)

## "FIXME" visually compare with numerical derivatives (e.g. from splines).
