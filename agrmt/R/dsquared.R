dsquared <- function(V) {
# calculate ordinal concentration, following Blair & Lacy 2000
# argument: V = frequency vector
k <- length(V) # number of categories
n <- sum(V)    # number of cases
V <- V/n       # standardizing
# from 1 to k-1(Fi -.5)^2
dsq <- sum(sapply(1:(k-1), function(x) (sum(V[1:x]) - .5)^2))
return(dsq)
}
