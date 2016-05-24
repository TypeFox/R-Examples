
# Infinitely differentiable transformations of R into a bounded or half-bounded interval

#Transformation from (-Inf, +Inf) to (-1, 0)
Trans.m10 <- function(x)
  -1/(1+exp(-x))
#Inverse
InvT.m10 <- function(x) 
  log(-x/(1+x))

#Transformation from (-Inf, +Inf) to (0, 1)
Trans.01 <- function(x)
  1/(1+exp(-x))
#Inverse
InvT.01 <- function(x)
  log(x/(1-x))

#Transformation from (-Inf, +Inf) to (0, +Inf)
Trans.0Inf <- function(x)
  exp(x)
#Inverse
InvT.0Inf <- function(x)
  log(x)

#Transformation from (-Inf, +Inf) to (1, +Inf)
Trans.1Inf <- function(x)
  1+exp(x)
#Inverse
InvT.1Inf <- function(x)
  log(x-1)

