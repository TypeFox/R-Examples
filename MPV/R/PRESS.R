PRESS <- function(x) {
  sum(resid(x)^2/(1-lm.influence(x)$hat)^2)
}
