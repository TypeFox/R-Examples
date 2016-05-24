########################################################################
## Stochastic Cellular Automaton
########################################################################
require(simecol)
CA <- new("gridModel",
  main = function(time, z, parms) {
    with(parms, {
      nb    <- eightneighbours(z)
      pgen  <- 1 - (1 - pbirth)^nb
      zgen  <- ifelse(z == 0 & runif(z) < pgen,         1,     0)
      zsurv <- ifelse(z >= 1 & runif(z) < (1 - pdeath), z + 1, 0)
      zgen + zsurv
    })
  },
  parms = list(pbirth = 0.02, pdeath = 0.01),
  times = c(from = 1, to = 50, by = 1),
  init = matrix(0, nrow = 80, ncol = 80),
  solver = "iteration",
  initfunc = function(obj) {
    init(obj)[38:42, 38:42] <- 5 # deterministic seed in the middle of the grid
    obj
  }
)

