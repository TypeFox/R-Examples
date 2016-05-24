################################################################################
##
## Example on using compiled code *without* simecol
##
## more information:
##   http://cran.r-project.org/web/packages/deSolve/vignettes/compiledCode.pdf
##
## thomas.petzoldt@tu-dresden.de
##
################################################################################



library("deSolve")

# compile C++ code within R
# (requires installed compiler)
# on Windows: http://www.murdoch-sutherland.com/Rtools/
system("R CMD SHLIB clotka.c")


modeldll <- dyn.load("clotka.dll")

parms <- c(k1 = 0.1, k2 = 0.1, k3=0.1) # parameters
init    <- c(prey = 1.0, pred = 0.5)   # initial conditions
times <- seq(0, 500, 0.1)              # output times


out <- lsoda(init, times, "dlotka", parms = parms, dllname = "clotka", nout = 2)
out <- as.data.frame(out)

plot(out$time, out$prey, type="l")
lines(out$time, out$pred, col="red")

dyn.unload(as.character(modeldll[2]))
