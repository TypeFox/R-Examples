sigmoid <- function(x) 1/(1 + exp(-x))

dsigmoid <- function(x) sigmoid(x) * (1 - sigmoid(x))

d2sigmoid <- function(x) dsigmoid(x) * (1 - 2 * sigmoid(x))
