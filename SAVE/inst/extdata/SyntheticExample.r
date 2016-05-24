#############
# Simulated example
#############

library("SAVE")

#######
# load data
#######

data(synthfield, package = "SAVE")
data(synthmodel, package = "SAVE")


##############
# create the SAVE object which describes the problem and
# compute the corresponding mle estimates
##############

synth <- SAVE(response.name = "y", controllable.names = "x", calibration.names = "v", 
              field.data = synthfield, model.data = synthmodel, mean.formula = ~1 + x, bestguess = list(v = 1.5))

##############
# Bayesian fit
##############

set.seed(0)
synth <- bayesfit(object = synth, prior = uniform(var.name = "v", lower = 0, upper = 3), 
                  n.iter = 20000)

summary(synth)

# Fig 9 Histogram of the posterior for v
plot(synth, option = "calibration")

# Fig 10 Histogram of the precisions
plot(synth, option = "precision")

#############
# validate at a grid of x points
##############

xnew <- data.frame(x = seq(from = 0.05, to = 3.05, length = 25))

valsynth <- validate(object = synth, newdesign = xnew, n.burnin = 100)

# summary and plot of the validation process
summary(valsynth)
plot(valsynth)

# Fig 11 A plot to compare the bias corrected prediction, the actual reality
# curve the computer model evaluated at the posterior mean and at the least
# square estimate
par(mfrow = c(1, 1))
a <- 0
b <- 5
av.real <- (valsynth@validate)[, "bias.corrected"]
delta <- (valsynth@validate)[, "tau.bc"]

# Bias corrected prediction and tolerance bars
plot(xnew$x, av.real, ty = "n", ylim = c(a, b), main = "Predictions", xlab = "x", 
     ylab = "y")
lines(xnew$x, av.real, lty = 1)
lines(xnew$x, av.real + delta, lty = 2)
lines(xnew$x, av.real - delta, lty = 2)
points(synthfield$x, synthfield$y, pch = "*")

# Representation of reality
yR <- function(x) {
  3.5 * exp(-1.7 * x) + 1.5
}

lines(xnew$x, yR(xnew$x), col = 2)

# Representation of Computer model at the least square estimate
yM <- function(x) {
  5 * exp(-x[1] * x[2])
}

vls <- 0.63
tmp <- cbind(xnew$x, vls)

lines(xnew$x, apply(tmp, 1, FUN = yM), col = 3)

# Representation of Computer model at the posterior mean
vmean <- mean(synth@mcmcsample[-(1:100), "v"])
tmp <- cbind(xnew$x, vmean)

lines(xnew$x, apply(tmp, 1, yM), col = 4)

