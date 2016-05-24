########################################################################
## Deterministic Cellular Automaton
##   Conway's Game of Life
########################################################################

set.seed(23)

conway <- new("gridModel",
  main = function(time, init, parms) {
    x      <- init
    nb     <- eightneighbours(x)
    surviv <- (x >  0 & (nb %in% parms$srv))
    gener  <- (x == 0 & (nb %in% parms$gen))
    x <- matrix((surviv + gener) > 0, nrow = nrow(init))
    return(x)
  },
  parms  = list(srv = c(2, 3), gen = 3),
  times  = 1:17,
  init   = matrix(round(runif(1000)), ncol = 40),
  solver = "iteration"
)

plot(sim(conway), axes = FALSE)


