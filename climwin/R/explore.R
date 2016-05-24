#'Visualise the weight distribution for given parameter values
#' 
#'Create a plot of the weibell or generalised extreme values (GEV) distribution
#'for given values of shape, scale and location parameters. Used to determine 
#'initial parameter values for \code{\link{weightwin}}.
#'@param shape A parameter that determines the shape of the distribution. 
#' Should be greater than 0.
#'@param scale A parameter that determines the scale of the distribution. 
#' Should be greater than 0.
#'@param loc A parameter that determines the location of the distribution. 
#' Should be less than or equal to 0.
#'@param weightfunc Choose whether to use a weibull ("W") or GEV ("G")
#' distribution.
#'@return ExploreWeight will return an example plot of the distribution using
#' given parameter values. This can be used to select the initial parameter 
#' values for \code{\link{weightwin}}
#'@author Martijn van de Pol and Liam D. Bailey
#'@examples
#'# Test a weibull distribution
#'
#'explore(shape = 3, scale = 0.2, loc = 0, weightfunc = "W")
#'
#'# Test a GEV distribution
#'
#'explore(shape = 3, scale = 5, loc = -5, weightfunc = "G")
#'
#'@export



explore <- function(shape = 1, scale = 1, loc = 0, weightfunc = "W"){
  par(mfrow = c(1, 1))
  duration <- 365
  j        <- seq(1:duration) / duration
  k        <- seq(-10, 10, by = (2 * 10 / duration))
  
  if (weightfunc == "W"){
    if (shape <= 0){
      stop("Weibull shape parameter should be >0")
    }
    
    if (scale <= 0){
      stop("Weibull scale parameter should be >0")
    }
    
    if (loc > 0){
      stop("Weibull location parameter should be <=0")
    }
    
    weight <- weibull3(x = j[1:duration], shape = shape, scale = scale, location = loc)
    weight[is.na(weight)] <- 0
    
    if (sum(weight) == 0){
      weight <- weight + 1
    }
    
    title <- c("Weibull", paste("shape=", shape, "scale=", scale, "location=", loc))
    plot((weight / sum(weight)), type = "l", ylab = "weight", xlab = "timesteps", main = title, xaxt = 'n')
  }
  
  if (weightfunc == "G"){
    if (scale <= 0){
      stop("GEV scale parameter should be >0")
    }
    weight <- dgev(k[1:duration], loc = loc, scale = scale, shape = shape, log = FALSE)
    weight[is.na(weight)] <- 0
    
    if (sum(weight) == 0){
      weight <- weight + 1
    }
    title <- c("GEV", paste("shape=", shape, "scale=", scale, "location=", loc))
    plot((weight / sum(weight)), type = "l", ylab = "weight", xlab = "timesteps", main = title, xaxt = 'n')
  }
}