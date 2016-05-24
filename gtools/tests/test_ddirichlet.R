## Regression test for bug reported by John Nolan:
##
## Whenever the pair (x[i], a[i]) == (0,1), NA would be returned, due
## to an internal computation of ( 0 * -Inf ) => NaN
##
## The code now checks for this particular issue and sets the value of
## ( 0 * -Inf ) to 0, which is correct for this calculation.
##
library(gtools)

x = c(0,0,1)
alpha = c(1,2,3)

stopifnot( ddirichlet(x=x, alpha=alpha) == 0 )

