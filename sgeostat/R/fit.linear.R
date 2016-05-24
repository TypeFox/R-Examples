assign("fit.linear",
function(v.object,type='c',plot.it=FALSE,iterations=1,c0=0,cl=1){
  # This program fits a univariate linear model to an empirical variogram
  # estimate.  The SEMI variogram model is the model fit...

  if (!inherits(v.object,"variogram")) stop('v.object must be of class, "variogram".\n')

# Interpret the "type" argument...
  if (type == 'c') empgamma <- v.object$classic /2
  else if (type == 'r') empgamma <- v.object$robust /2
  else if (type == 'm') empgamma <- v.object$med /2
  else stop("type must be 'classic', 'robust', or 'median'.\n")

# Set up the variogram function...
  linear.v <- function(h,parameters)
                   ifelse(h==0,0,
                   parameters[1] + parameters[2]*h)

# Get the number of observations and the bins (h's)...
  numobs <- v.object$n
  h      <- v.object$bins

# If any numobs is 0, get rid of that lag...
  empgamma <- empgamma[numobs>0]
  h <- h[numobs>0]
  numobs <- numobs[numobs>0]

  if(iterations>0){
    fit <- lsfit(h,empgamma)
    parameters <- fit$coef
  }
  else
    parameters <- c(c0,cl)
  
  names(parameters)<-c("nugget","slope")

  v.m.object <- list(parameters=parameters,
                     model= linear.v
                    )
  attr(v.m.object,'class') <- 'variogram.model'
  attr(v.m.object,'type') <- 'linear'
  if (plot.it)
    plot(v.object,var.mod.obj=v.m.object)

  return(v.m.object)
  
})
