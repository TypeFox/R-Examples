eapEst <-
function( resp,                             # The vector of responses
          params,                           # The item parameters
          range = c(-6, 6),                 # The limits of integration
          mod = c("brm", "grm"),            # The model                        
          ddist = dnorm, quad = 33, ... )   # The prior distribution stuff:
{
  
# First turn params into a matrix:
  params <- rbind(params)

# And turn response into a matrix:
  resp <- { if( dim(params)[1] > 1 ) rbind(resp)   # ... --> turn it into a multi-column matrix,
            else                     cbind(resp) } # ... --> or a 1-column matrix

  class(params) <- c(mod, "matrix")

#~~~~~~~~~~~~~~~~~#
# Argument Checks #
#~~~~~~~~~~~~~~~~~#

# Make sure that the arguments are OK:

## 1 ## (Make sure that params and resp are ALL numeric)
  if( mode(params) != "numeric" )
    stop( "params need to be numeric" )
    
  if( !is.null(resp) & mode(resp) != "numeric" )
    stop( "resp needs to be numeric" )

## 2 ## (Make sure that the dimensions of params and response are equal)
  if( !is.null(resp) & ( dim(resp)[2] != dim(params)[1] ) )
    stop( "number of params does not match the length of resp" )


#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Expected A Posteriori Est #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Indicate the lower/upper boundary of integration:
  if( is.null(range) )
    range <- c(-6, 6)

  l <- range[1]; u <- range[2]
  
# Make sure that the lower/upper points of integration is within the legal range:
  tmp.x <- seq(l, u, by = .01)
  tmp.y <- ddist(x = tmp.x, ... )
  
# Usually, if the density doesn't exist, it prints "0", but it might print NA:
  l <- min( tmp.x[signif(tmp.y) != 0 & !is.na(tmp.y)] )
  u <- max( tmp.x[signif(tmp.y) != 0 & !is.na(tmp.y)] )
  
  est <- NULL # a vector for estimates
  sem <- NULL # a vector for posterior var/sds

# The Likelihood function used in the EAP integration:
  LikFun <- function( ... ) exp( get(paste("logLik.", mod, sep = "") )( ... ) )

# For each person:
  for( i in 1:dim(resp)[1] ){
  	
  	resp_i <- resp[i, ]
    
# Set up the numerators and denominator of integral:
    eap.den <- function( theta )           LikFun(u = resp_i, theta = theta, params = params) * ddist(x = theta, ... )
    eap.num <- function( theta ) theta   * LikFun(u = resp_i, theta = theta, params = params) * ddist(x = theta, ... )
    sem.num <- function( theta ) theta^2 * LikFun(u = resp_i, theta = theta, params = params) * ddist(x = theta, ... )
    
# Set up integration function:
    integrate.eap <- function(f) integrate.s(f, min = l, max = u, quad = quad, method = "S2")
    
# And integrate (sem is the variance at this point):
     est[i] <- ( integrate.eap(eap.num) / ( const <- integrate.eap(eap.den) ) )
     sem[i] <- ( integrate.eap(sem.num) / const ) - est[i]^2
  
  } # END FOR LOOP
  
# Round the estimated value to three/four? decimal places and find the sem:          
  est <- round(est, digits = 4)
  sem <- sqrt(sem)
  
# And pull out the information:
  info <- get(paste("FI.", mod, sep = ""))( params = params,
                                            theta = est,
                                            type = "observed",
                                            resp = resp )$test
  
  return( list( theta = est, info = info, sem = sem ) )
  
} # END eapEst FUNCTION