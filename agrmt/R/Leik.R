Leik <- function (V) {
   # Leik (1966) A Measure of Ordinal Consensus
   # Argument: V = vector containing frequnecy distribution ("frequency vector")
   n <- sum(V)     # sum of all values
   m <- length(V)  # number of categories
   P <- V/n        # percentages
   R <- cumsum(P)  # cumulative frequency distribution
   SDi <- sapply(1:m, function(x) ifelse(R[x] <= .5, R[x], 1-R[x]))  # Differences, eqn.1
   maxSDi <- .5*(m-1)             # Maximum possible SDi given n, eqn.2
   D <- sum(SDi)/maxSDi           # ordinal dispersion; standardized, eqn.3
   return(D)
   }
