assign("fit.gaussian",
function(v.object,c0=NULL,cg=NULL,ag=NULL,type='c',iterations=10,
         tolerance=.000001,echo=FALSE,plot.it=FALSE,weighted=TRUE){
  # This program fits a univariate gaussian model to an empirical variogram
  # estimate.  The SEMI variogram model is the model fit...

  if (!inherits(v.object,"variogram")) stop('v.object must be of class, "variogram".\n')

  if (is.null(c0) & is.null(cg) & is.null(ag)) estimate.initial <- TRUE
  else estimate.initial <- FALSE

  if ( !estimate.initial & (is.null(c0) | is.null(cg) | is.null(ag)))
    stop('c0, cg, and ag must all be entered.\n')

# Interpret the "type" argument...
  if (type == 'c') empgamma <- v.object$classic /2
  else if (type == 'r') empgamma <- v.object$robust /2
  else if (type == 'm') empgamma <- v.object$med /2
  else stop("type must be 'c', 'r', or 'm'.\n")

# Set up the first derivative functions for each of the parameters...
# dc0 = 1
  dcg <- function(h,ag) { return(1-exp(-h^2/ag^2)) }
  dag <- function(h,cg,ag) { return(-2*cg*h^2/(ag^3)*
                              exp(-h^2/ag^2))}

# Set up the variogram function...
  gaussian.v <- function(h,parameters)
                   ifelse(h==0,0,
                   parameters[1] + parameters[2]*(1-exp(-h^2/parameters[3]^2)))

# Get the number of observations and the bins (h's)...
  numobs <- v.object$n
  h      <- v.object$bins

# If any numobs is 0, get rid of that lag...
  empgamma <- empgamma[numobs>0]
  h <- h[numobs>0]
  numobs <- numobs[numobs>0]

# Start the sums of squares at 0...
  rse <- 0

# Begin iterations...
  parameters <- c(c0,cg,ag)
  if(iterations>0){
    cat('Initial parameter estimates: ',parameters,'\n')
    loop <- TRUE
    converge <- FALSE
    i <- 1
  } else {
    loop <- FALSE
    converge <- FALSE
  }

# Plot it before we start if requested...
  if (plot.it) {
    v.m.object <- list(parameters=parameters,
                       model= gaussian.v
                      )
    attr(v.m.object,'class') <- 'variogram.model'
    attr(v.m.object,'type') <- 'gaussian'

    plot(v.object,var.mod.obj=v.m.object,type=type)
  }

  while (loop) {
    cat('Iteration:',i,'\n')
# establish the Y vector...
    y <- (empgamma - gaussian.v(h,parameters))

# establish the x matrix...
    xmat <- cbind(rep(1,length(h)),
                  dcg(h,parameters[3]), 
                  dag(h,parameters[2],parameters[3]))

# establish the weights (Cressie, p. 99)...
    if(weighted) {
      w <- numobs/(gaussian.v(h,parameters))^2
      }
    else {
#      w <- (1:length(numobs))
      w <- rep(1,length(numobs))
      }
   
    if(echo) cat('  X matrix:\n')
    if(echo) print(cbind(y,xmat,w))
    if(echo) cat('\n\n')    

    fit <- lsfit(xmat,y,wt=w,intercept=FALSE)

# calculate the new parameter estimates...
    parameters.old <- parameters
    parameters <- fit$coef + parameters
    parameters <- ifelse(parameters>0,parameters,.000001)
    cat('Gradient vector: ',fit$coef,'\n')
    cat('New parameter estimates: ',parameters,'\n\n')

# Check for convergence, see if the sum of squares has converged
    rse.old <- rse
    rse <- sum(fit$residuals^2)
    rse.dif <- rse-rse.old

# Check for convergence of parmeters...
    parm.dist <- sqrt(sum((parameters-parameters.old)^2))

    cat('rse.dif = ',rse.dif,'(rse =',rse,')  ;  parm.dist = ',parm.dist,'\n\n')
#    cat('rse.dif = ',rse.dif,'(rse =',rse,')\n\n')
#    if(rse.dif < tolerance & parm.dist < tolerance) {
    if(abs(rse.dif) < tolerance) {
      loop <- FALSE
      converge <- TRUE
      cat('Convergence achieved by sums of squares.\n')
    }

    i <- i+1
    if (i>iterations) {
      loop <- FALSE
    }
    v.m.object <- list(parameters=parameters,
                       model= gaussian.v
                      )
    attr(v.m.object,'class') <- 'variogram.model'
    attr(v.m.object,'type') <- 'gaussian'
    if (plot.it)
      plot(v.object,var.mod.obj=v.m.object,type=type)
  }
  if (converge) 
    cat('Final parameter estimates: ',parameters,'\n\n')
  else
    cat('Convergence not achieved!\n')
  
  names(parameters)<-c("nugget","sill","range")
  v.m.object <- list(parameters=parameters,
                     model= gaussian.v
                    )
  attr(v.m.object,'class') <- 'variogram.model'
  attr(v.m.object,'type') <- 'gaussian'

  return(v.m.object)
  
})
