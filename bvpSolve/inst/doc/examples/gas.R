## =============================================================================
## Emden's equation
## Spherical body of gas
## y''+2/x y' +y^n=0
## y'0=0, y(1)=sqrt(3/4)
##
## becomes:
## dy=y2
## dy2=-2/x y2 -y^n
## =============================================================================

require(bvpSolve)

Gas<-function(x, y, par)
{
  dy1 <- y[2]
  dy2 <- -2/x*y[2]-(y[1]^n)
  if (x < 1e-10) dy2<-0     # quick and dirty
  list(c(dy1, dy2))
}

x  <-seq(0, 1, by = 0.01)
n  <- 5

## =============================================================================
## 1. shooting method
## =============================================================================
print(system.time(
sol <- bvpshoot(func = Gas, yini = c(y = NA, dy = 0),
    yend = c(sqrt(3/4), NA), x = x, guess = 0)
))

plot(sol)
# add analytical solution
curve(1/sqrt(1+(x^2)/3), type = "l", add = TRUE)

## =============================================================================
## 2. bvptwp method
## =============================================================================

print(system.time(
Sol <- bvptwp(func = Gas, yini = c(NA, 0), yend = c(sqrt(3/4), NA),
             x = x, parms = NULL)
))
lines(Sol, col = "red")

