##
## runit-gp.r - unit tests of the (low level) genetic programming function
##

## test the basic functionality of the GP loop
test.basic_gp <- function() {
  sinusfitness <- makeFunctionFitnessFunction(sin, -pi, pi, steps = 256, indsizelimit = 32)
  pop1 <- geneticProgramming(sinusfitness, stopCondition=makeStepsStopCondition(100))
  checkTrue(length(pop1$population) == 500)
}
