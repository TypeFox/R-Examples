################################################
## FRACTAL chaotic systems 
##
##   henon
##   lorenz.ode
##
################################################

###
# henon
###

"henon" <- function(start=rnorm(2), a=1.4, b=0.3, n.sample=2000, n.transient=10)
{
  # check argument types and length
  checkVectorType(start,"numeric")
  if (length(start) != 2)
    stop("start vector must contain two numeric values")
  checkScalarType(a,"numeric")
  checkScalarType(b,"numeric")
  checkScalarType(n.sample,"integer")
  checkScalarType(n.transient,"integer")
  
  # check argument values
  if (n.sample < 1)
    stop("n.sample must be positive")
  if (n.transient < 1)
    stop("n.transient must be positive")

  # initialize variables
  n.sample <- n.sample + n.transient
  y <- x <- vector(mode="numeric", length=n.sample)
  x[1] <- start[1]
  y[1] <- start[2]

  # iterate map
  for (n in seq(2, n.sample)){

    x[n] <- a - x[n - 1]^2 + b * y[n - 1]
    y[n] <- x[n - 1]
  }

  # remove transient
  transient <- seq(n.transient)

  # return output
  list(x=x[-transient], y=y[-transient])
}

###
# lorenz.ode
###

"lorenz.ode" <- function(x, sigma=10, r=28, b=8/3)
{
  # initialize variables
  if (is.null(x))
    x <- rnorm(3)

  # check argument types and length
  checkVectorType(x,"numeric")
  if (length(x) != 3)
    stop("x vector must contain two numeric values")
  checkScalarType(sigma,"numeric")
  checkScalarType(r,"numeric")
  checkScalarType(b,"numeric")

  # return the next state of the Lorenz system
  c(sigma * (x[2] - x[1]), x[1] * (r - x[3]) - x[2], -b * x[3] + x[1] * x[2])
}
