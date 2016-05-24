options(digits=5)
library(gee)
load("data_for_gee_binomial.RData")
try(gee(yb ~ x, family = binomial, id = id, R = R, cor = "fixed"))
## infinite looped in gee 4.3-14


