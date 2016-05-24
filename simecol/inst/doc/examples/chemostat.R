chemostat <- new("odeModel",
  main = function(time, init, parms, inputs = NULL) {
    with(as.list(c(init, parms)), {
      mu  <- vm * S / (km + S)            # Monod equation
      dx1 <- mu * X - D * X               # cells, e.g. algae
      dx2 <-  D *(S0 - S) - 1/Y * mu * X  # substrate, e.g. phosphorus
      list(c(dx1, dx2))
    })
  },
  parms = c(
    vm = 1.0,           # max growth rate, 1/d
    km = 2.0,           # half saturation constant, mumol / L
    Y  = 100,           # cells /mumol Substrate
    D  = 0.5,           # dilution rate, 1/d
    S0 = 10             # substrate in inflow, mumol / L
  ),
  times = c(from = 0, to = 40, by = .5),
  init  = c(X = 10, S = 10), # cells / L; Substrate umol / L
  solver = "lsoda"
)
