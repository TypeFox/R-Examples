IKL <- function( params,      # parameters over which to calculate
                 theta,       # a single estimate of theta
                 delta = .1,  # the indifference region specification
                 quad = 33 )  # the number of quadrature points
{
  
# This function finds Integrated KL-Information:

#~~~~~~~~~~~~~~~~~~~~~~~#
# Bounds of Integration #
#~~~~~~~~~~~~~~~~~~~~~~~#

# Indicate the lower/upper bounds of integration and density:
  range <- theta + c(-1, 1)*delta
    
  l <- range[1]; u <- range[2]
  
# Set up the bounds of integration:
  X       <- seq(l, u, length = quad)  # region of integration
  X.theta <- (X + theta)/2             # shifting center of indifference region
  X.delta <- X - X.theta               # half-width of indifference region
  
# The integral is int_{theta' - delta}^{theta' + delta}KL(theta | theta')
#  --> So we have several KL information thingies
#  --> The center of each KL information: theta'' = (theta + theta')/2
#  --> The half-width of the corresponding indifference region: abs(theta'' - theta)

# Note: 1) Shifted KL information as we go from theta' - delta to theta' + delta
#       2) X.theta - X.delta will always/must always equal theta for all points.

#~~~~~~~~~~~~~~~~#
# KL-Information #
#~~~~~~~~~~~~~~~~#
               
# The KL-Information for all of the items thus far across range:
  Y <- KL(params = params, theta = X.theta, delta = X.delta)$item

#~~~~~~~~~~~~~~~#
# And Integrate #
#~~~~~~~~~~~~~~~#

  info <- apply(Y, MARGIN = 2, FUN = integrate.q, x = X)
  
  return( list(item = info) )
  
} # END IKL FUNCTION
