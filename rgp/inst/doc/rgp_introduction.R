### R code from vignette source 'rgp_introduction.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: installation (eval = FALSE)
###################################################
## install.packages("rgp")


###################################################
### code chunk number 2: library
###################################################
library("rgp")


###################################################
### code chunk number 3: sine_unction_set_1
###################################################
functionSet1 <- functionSet("+", "*", "-")


###################################################
### code chunk number 4: sine_input_variable_set_1
###################################################
inputVariableSet1 <- inputVariableSet("x")


###################################################
### code chunk number 5: sine_constant_factory_set_1
###################################################
constantFactorySet1 <- constantFactorySet(function() rnorm(1))


###################################################
### code chunk number 6: sine_fitness_function_1
###################################################
interval1 <- seq(from = -pi, to = pi, by = 0.1)
fitnessFunction1 <- function(f) rmse(f(interval1), sin(interval1))


###################################################
### code chunk number 7: sine_gp_run_1 (eval = FALSE)
###################################################
## set.seed(1)
## gpResult1 <- geneticProgramming(functionSet = functionSet1,
##                                 inputVariables = inputVariableSet1,
##                                 constantSet = constantFactorySet1,
##                                 fitnessFunction = fitnessFunction1,    
##                                 stopCondition = makeTimeStopCondition(5 * 60))


###################################################
### code chunk number 8: sine_best_solution_1 (eval = FALSE)
###################################################
## bestSolution1 <- gpResult1$population[[which.min(gpResult1$fitnessValues)]]


###################################################
### code chunk number 9: sine_best_solution_1_plot (eval = FALSE)
###################################################
## plot(y = bestSolution1(interval1), x = interval1, type = "l",
##   lty = 1, xlab = "x", ylab = "y")
## lines(y = sin(interval1), x = interval1, lty = 2)


###################################################
### code chunk number 10: pendulum
###################################################
makeDampedPendulum <- function(A0 = 1, g = 9.81, l = 0.1, phi = pi, gamma = 0.5) {
  omega <- sqrt(g/l)
  function(t) A0 * exp(-gamma * t) * cos(omega * t + phi)
}


###################################################
### code chunk number 11: pendulum_examples
###################################################
pendulum1 <- makeDampedPendulum(l = 0.5)
pendulum2 <- makeDampedPendulum(l = 1.2, A0 = 0.5)


###################################################
### code chunk number 12: pendulum_example_plots
###################################################
interval1 <- seq(from = 0, to = 10, by = 0.05)
plot(y = pendulum1(interval1), x = interval1, type = "l",
  lty = 1, xlab = "t", ylab = "deflection")
lines(y = pendulum2(interval1), x = interval1, lty = 2)


###################################################
### code chunk number 13: rgp_introduction.Rnw:483-484
###################################################
interval1 <- seq(from = 0, to = 10, by = 0.05)
plot(y = pendulum1(interval1), x = interval1, type = "l",
  lty = 1, xlab = "t", ylab = "deflection")
lines(y = pendulum2(interval1), x = interval1, lty = 2)


###################################################
### code chunk number 14: pendulum_data
###################################################
  xs1 <- seq(from = 1, to = 10, length.out = 512)
  pendulum1Data <- data.frame(time = xs1,
    deflection = pendulum1(xs1) + rnorm(length(xs1), sd = 0.01))


###################################################
### code chunk number 15: pendulum_gp_run (eval = FALSE)
###################################################
## modelSet1 <- symbolicRegression(deflection ~ time, data = pendulum1Data,
##                                 stopCondition = makeTimeStopCondition(2 * 60))


###################################################
### code chunk number 16: pendulum_result_analysis (eval = FALSE)
###################################################
## bestModel1 <- modelSet1$population[[which.min(modelSet1$fitnessValues)]]
## plot(y = bestModel1(xs1), x = xs1, type = "l",
##   lty = 1, xlab = "x", ylab = "y")
## lines(y = pendulum1(xs1), x = xs1, lty = 2)


###################################################
### code chunk number 17: typed_parity
###################################################
parity <- function(x) {
  numberOfOnes <- sum(sapply(x, function(bit) if (bit) 1 else 0))
  numberOfOnes %% 2 != 0
}


###################################################
### code chunk number 18: typed_parity3
###################################################
parity3 <- function(x1, x2, x3) parity(c(x1, x2, x3))


###################################################
### code chunk number 19: typed_parity_fitness
###################################################
parityFitnessFunction <- makeBooleanFitnessFunction(parity3)


###################################################
### code chunk number 20: typed_boolean_constant_set
###################################################
booleanConstantFactory <- function() runif(1) > .5
booleanConstantSet <- constantFactorySet(
  "booleanConstantFactory" %::% (list() %->% st("logical")))


###################################################
### code chunk number 21: typed_boolean_function_set
###################################################
booleanFunctionSet <- functionSet(
  "&" %::% (list(st("logical"), st("logical")) %->% st("logical")),
  "|" %::% (list(st("logical"), st("logical")) %->% st("logical")),
  "!" %::% (list(st("logical")) %->% st("logical")))


###################################################
### code chunk number 22: typed_boolean_input_variable_set
###################################################
booleanInputVariableSet <- inputVariableSet(
  "x1" %::% st("logical"), 
  "x2" %::% st("logical"), 
  "x3" %::% st("logical"))


###################################################
### code chunk number 23: typed_gp_run (eval = FALSE)
###################################################
## typedGpResult1 <- typedGeneticProgramming(parityFitnessFunction, st("logical"),
##   functionSet = booleanFunctionSet,
##   inputVariables = booleanInputVariableSet,
##   constantSet = booleanConstantSet,
##   stopCondition = makeTimeStopCondition(30))


###################################################
### code chunk number 24: typed_gp_result_analysis (eval = FALSE)
###################################################
## bestFunction1 <- typedGpResult1$population[[which.min(typedGpResult1$fitnessValues)]]


###################################################
### code chunk number 25: spot_installation (eval = FALSE)
###################################################
## install.packages("SPOT")


###################################################
### code chunk number 26: spot_rgp_ (eval = FALSE)
###################################################
## x1 <- seq(0, 4 * pi, length.out = 201)
## x2 <- seq(0, 4 * pi, length.out = 201)
## y <- sin(x1) + cos(2 * x2)
## data1 <- data.frame(y = y, x1 = x1, x2 = x2)
## 
## result1 <- symbolicRegression(y ~ x1 + x2,
##   data = data1,
##   populationSize = populationSize,
##   selectionFunction = makeTournamentSelection(tournamentSize = tournamentSize),
##   functionSet = arithmeticFunctionSet,
##   stopCondition = makeTimeStopCondition(time))
## 
## bestFitness <- min(sapply(result1$population, result$fitnessFunction))


###################################################
### code chunk number 27: spot_run (eval = FALSE)
###################################################
## library("SPOT")
## confPath <- find.package("SPOT")
## confPath <- file.path(confPath, "demo17Rgp")
## confFile <- file.path(confPath, "rgp0001.conf")
## spotConfig <- spot(confFile)


