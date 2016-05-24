\dontrun{

##############################################################################################
## Estimate parameter values from experimental data
##############################################################################################

library(deSolve)

# O2 + O <-> O3
# diff = O2 - O3
# build_O3 = const.
f <- c(
  O3 = " build_O3 * O2 * O - decay_O3 * O3",
  O2 = "-build_O3 * O2 * O + decay_O3 * O3",
  O  = "-build_O3 * O2 * O + decay_O3 * O3"
)

# Compute sensitivity equations and get attributes
f_s <- sensitivitiesSymb(f)
chi <- attr(f_s, "chi")
grad <- attr(f_s, "grad")
forcings <- attr(f_s, "forcings")

# Generate ODE function
func <- funC(f = c(f, f_s, chi, grad), forcings = forcings, fcontrol = "nospline")

# Initialize times, states, parameters
times <- seq(0, 15, by = .1)
yini <- c(O3 = 0, O2 = 2, O = 2.5)
yini_s <- attr(f_s, "yini")
yini_chi <- c(chi = 0)
yini_grad <- rep(0, length(grad)); names(yini_grad) <- names(grad)
pars <- c(build_O3 = .2, decay_O3 = .1)

# Initialize forcings (the data)
data(oxygenData)
forcData <- data.frame(time = oxygenData[,1],
                       name = rep(colnames(oxygenData[,-1]), each=dim(oxygenData)[1]),
                       value = as.vector(oxygenData[,-1]))
forc <- setForcings(func, forcData)

# Solve ODE
out <- odeC(y = c(yini, yini_s, yini_chi, yini_grad), 
            times = times, func = func, parms = pars, forcings = forc,
            method = "lsodes")

# Plot solution
par(mfcol=c(1,2))
t <- out[,1]
M1 <- out[,2:4]
M2 <- out[,names(grad)]
tD <- oxygenData[,1]
M1D <- oxygenData[,2:4]

matplot(t, M1, type="l", lty=1, col=1:3, xlab="time", ylab="value", main="states")
matplot(tD, M1D, type="b", lty=2, col=1:3, pch=4, add=TRUE)
legend("topright", legend = names(f), lty=1, col=1:3)
matplot(t, M2, type="l", lty=1, col=1:5, xlab="time", ylab="value", main="gradient")
legend("topleft", legend = names(grad), lty=1, col=1:5)

# Define objective function
obj <- function(p) {
  out <- odeC(y = c(p[names(f)], yini_s, yini_chi, yini_grad), 
              times = times, func = func, parms = p[names(pars)], 
	      forcings = forc, method="lsodes")
  
  value <- as.vector(tail(out, 1)[,"chi"])
  gradient <- as.vector(tail(out, 1)[,paste("chi", names(p), sep=".")])
  hessian <- gradient%*%t(gradient)
  
  return(list(value = value, gradient = gradient, hessian = hessian))
}

# Fit the data
myfit <- optim(par = c(yini, pars), 
               fn = function(p) obj(p)$value, 
               gr = function(p) obj(p)$gradient,
               method = "L-BFGS-B", 
               lower=0,
               upper=5)

# Model prediction for fit parameters
prediction <- odeC(y = c(myfit$par[1:3], yini_s, yini_chi, yini_grad), 
                   times = times, func = func, parms = myfit$par[4:5], 
		   forcings = forc, method = "lsodes")

# Plot solution
par(mfcol=c(1,2))
t <- prediction[,1]
M1 <- prediction[,2:4]
M2 <- prediction[,names(grad)]
tD <- oxygenData[,1]
M1D <- oxygenData[,2:4]

matplot(t, M1, type="l", lty=1, col=1:3, xlab="time", ylab="value", main="states")
matplot(tD, M1D, type="b", lty=2, col=1:3, pch=4, add=TRUE)
legend("topright", legend = names(f), lty=1, col=1:3)
matplot(t, M2, type="l", lty=1, col=1:5, xlab="time", ylab="value", main="gradient")
legend("topleft", legend = names(grad), lty=1, col=1:5)

}
