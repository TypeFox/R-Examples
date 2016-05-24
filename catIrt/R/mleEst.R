mleEst <-
function( resp,                         # The vector of responses
          params,                       # The item parameters
          range = c(-6, 6),             # The integer to maximize over
          mod = c("brm", "grm"),        # The model
          ... )
{

# First turn params into a matrix:
  params <- rbind(params)
   
# And turn response into a matrix:
  resp <- { if( dim(params)[1] > 1 ) rbind(resp)   # ... --> turn it into a multi-column matrix,
            else                     cbind(resp) } # ... --> or a 1-column matrix
  
#~~~~~~~~~~~~~~~~~#
# Argument Checks #
#~~~~~~~~~~~~~~~~~#

# Make sure that the arguments are OK:

## 1 ## (Make sure that params and resp are ALL numeric)
  if( mode(params) != "numeric" )
    stop("params need to be numeric")
    
  if( !is.null(resp) & mode(resp) != "numeric" )
    stop("resp needs to be numeric")

## 2 ## (Make sure that the dimensions of params and response are equal)
  if( !is.null(resp) & ( dim(resp)[2] != dim(params)[1] ) )
    stop("number of params does not match the length of resp")


#~~~~~~~~~~~~~~~~~~~~~~~~#
# Maximum Likelihood Est #
#~~~~~~~~~~~~~~~~~~~~~~~~#

# Indicate the lower/upper boundary of the search:
  if( is.null(range) )
    range <- c(-6, 6)

  l <- range[1]; u <- range[2]
  
  est <- NULL # a vector for estimates
  
# Then, maximize the loglikelihood function over that interval for each person:
  for( i in 1:dim(resp)[1] ){
    likFun <- paste("logLik.", mod, sep = "")
    est[i] <- optimize( get(likFun), lower = l, upper = u, maximum = TRUE,
                        u = resp[i, ], params = params,
                        type = "MLE" )$max
  } # END for LOOP
                
# Round the estimated value to three/four? decimal places:                      
  est <- round(est, digits = 4)
  
# And pull out the information as well as the SEM:
  info <- get(paste("FI.", mod, sep = ""))( params = params,
                                            theta = est,
                                            type = "observed",
                                            resp = resp )
  
  list( theta = est, info = info$test, sem = info$sem )
  
} # END mleEst FUNCTION