########################################################################
## Logistic growth (in ODE formulation)
########################################################################

logist <- new("odeModel",
  main = function (time, init, parms) {
    x <- init
    p <- parms
    dx1 <-   p["r"] * x[1] * (1 - x[1] / p["K"])
    list(c(dx1))
  },
  parms  = c(r=0.1, K=10),
  times  = seq(0, 100, 1),
  init   = c(population=0.1),
  solver = "lsoda"
)


