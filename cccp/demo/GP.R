##
## Demo for geometric program
## Example taken from:
## Beightler, C.S. and D.T. Phillips, Applied Geometric Programming,
## John Wiley and Sons, New York, NY, 1976.
##
## GP formulation:
## F0 = 0.44 * x_1^3 * x_2^-2 + 10 * x_1^-1 + 0.592 * x_1 * x_2^-3
## F1 = 8.62 * x_^1^-1 * x_2^3 <= 1.0 ; x_1, x_2 > 0
##
## Creating Problem
F0 <- matrix(c(3, -2, -1, 0, 1, -3), nrow = 3, ncol = 2, byrow = TRUE)
g0 <- log(c(0.44, 10, 0.592))
F1 <- matrix(c(-1, 3), nrow = 1, ncol = 2, byrow = TRUE)
g1 <- log(8.62)
ans <- gp(F0, g0, FList = list(F1), gList = list(g1))
ans
