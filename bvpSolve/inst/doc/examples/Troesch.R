## =============================================================================
## Troesch's equation
## =============================================================================
require(bvpSolve)

Troesch <- function (t, y, pars) {
  list( c(y[2],  mu*sinh(mu*y[1])) )
}

mu <- 5

# number of mesh points: 12
x <- seq(0, 1, len = 12)

Sol <- bvptwp(yini = c(y = 0, dy = NA), yend = c(1, NA),
          x = x, fun = Troesch)

plot(Sol)
