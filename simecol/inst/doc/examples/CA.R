########################################################################
## An Elementary Stochastic Cellular Automaton
########################################################################

CA <- new("gridModel",
    main = function(time, init, parms) {
      z     <- init
      nb    <- eightneighbors(z)
      pgen  <- 1 - (1 - parms$pbirth)^nb
      zgen  <- ifelse(z == 0 & 
                 runif(z) < pgen, 1, 0)
      zsurv <- ifelse(z >= 1 & 
                 runif(z) < (1 - parms$pdeath), 
                 z + 1, 0)
      zgen + zsurv
    },
    parms = list(pbirth = 0.02, pdeath = 0.01),
    times = c(from = 1, to = 50, by = 1),
    init = matrix(0, nrow = 40, ncol = 40),
    solver = "iteration"
)
init(CA)[18:22,18:22] <- 1

