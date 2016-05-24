\dontrun{

##############################################################################################
## Sensitivity analysis of ozone formation
##############################################################################################

library(deSolve)


# 0 -> 2A -> B -> C -> 0
f <- c(A = "-2*k1*A*A +1*buildA",
       B = "-1*k2*B +1*k1*A*A", 
       C = "-1*k3*C +1*k2*B")

forcings <- "buildA"
forcData <- data.frame(name = "buildA", time = c(0, 5, 10, 15), value = c(0, 1, 2, 0))

# Compute sensitivity equations
f_s <- sensitivitiesSymb(f, reduce = TRUE, inputs = forcings)
outputs <- c(attr(f_s, "outputs"), sum = "A + B + C")

# Generate ODE function
func <- funC(c(f, f_s), forcings = forcings, outputs = outputs, fcontrol = "nospline")

# Initialize times, states, parameters and forcings
times <- seq(0, 15, by = .1)
yini <- c(A = 0, B = 0, C = 0, attr(f_s, "yini"))
pars <- c(buildA = 1, k1 = .1, k2 = .01, k3 = .001)
forc <- setForcings(func, forcData)

# Solve ODE
out <- odeC(y = yini, times = times, func = func, parms = pars, forcings = forc)
class(out) <- "deSolve"
plot(out)
