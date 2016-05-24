\dontrun{

##############################################################################################
## Solve an optimal control problem:
##############################################################################################

library(bvpSolve)

# O2 + O <-> O3
# O3 is removed by a variable rate u(t)
f <- c(
  O3 = " build_O3 * O2 * O - decay_O3 * O3 - u * O3",
  O2 = "-build_O3 * O2 * O + decay_O3 * O3",
  O  = "-build_O3 * O2 * O + decay_O3 * O3"
)

# Compute adjoints equations and replace u by optimal input
f_a <- adjointSymb(f, states = c("O3"), inputs = "u")
inputs <- attr(f_a, "inputs")
f_tot <- replaceSymbols("u", inputs, c(f, f_a))
forcings <- attr(f_a, "forcings")

# Initialize times, states, parameters
times <- seq(0, 15, by = .1)
boundary <- data.frame(
  name = c("O3", "O2", "O", "adjO3", "adjO2", "adjO"),
  yini = c(0.5, 2, 2.5, NA, NA, NA),
  yend = c(NA, NA, NA, 0, 0, 0))

pars <- c(build_O3 = .2, decay_O3 = .1, eps = 1)

# Generate ODE function
func <- funC(f = f_tot, forcings = forcings, jacobian = "full", boundary = boundary)

# Initialize forcings (the objective)
forcData <- data.frame(time = times,
                       name = rep(forcings, each=length(times)),
                       value = rep(c(0.5, 0, 1, 1), each=length(times)))
forc <- setForcings(func, forcData)

# Solve BVP
out <- bvptwpC(x = times, func = func, parms = pars, forcings = forc)

# Plot solution
par(mfcol=c(1,2))
t <- out[,1]
M1 <- out[,2:4]
M2 <- with(list(uD = 0, O3 = out[,2], adjO3 = out[,5], eps = 1, weightuD = 1), 
           eval(parse(text=inputs)))

matplot(t, M1, type="l", lty=1, col=1:3, xlab="time", ylab="value", main="states")
abline(h = .5, lty=2)
legend("topright", legend = names(f), lty=1, col=1:3)
matplot(t, M2, type="l", lty=1, col=1, xlab="time", ylab="value", main="input u")
abline(h = 0, lty=2)

}
