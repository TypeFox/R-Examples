### inv mills ratio.
### Should accept all kind of values, including far off range.
lambda <- function(x) {
   as.vector(ifelse(x > -30, dnorm(x)/pnorm(x), -x))
                           # can we prove it?
}
