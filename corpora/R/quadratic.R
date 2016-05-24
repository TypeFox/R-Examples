##
## internal helper functions for quadratic equations
##

## solve quadratic equation a x^2 + b x + c = 0 for vectors of coefficients
solve.quadratic <- function (a, b, c, nan.upper=NA, nan.lower=NA) {
  d <- b * b - 4 * a * c                # discriminant
  data.frame(lower = ifelse(d < 0, rep(nan.upper, length(d)), (-b - sqrt(d)) / (2*a)),
             upper = ifelse(d < 0, rep(nan.lower, length(d)), (-b + sqrt(d)) / (2*a)))
}
