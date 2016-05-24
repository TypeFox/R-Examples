simIrt <-
function( theta = seq(-3, 3, by = .1), # a scalar/vector of theta values
          params,                      # the item parameters
          mod = c("brm", "grm") )      # the model (binary response model or graded response model for now)
{

#~~~~~~~~~~~~~~~~~#
# Argument Checks #
#~~~~~~~~~~~~~~~~~#

# We need to make sure that the arguments are OK:

## 1 ## (Make sure that theta is a numeric vector)
  if( !is.null( dim(theta) ) | mode(theta) != "numeric" )
    stop( "theta must be a numeric vector of person/simulee parameters" )
    
## 2 ## (Make sure that params is a numeric matrix)
  if( mode(params) != "numeric" & !inherits(params, "matrix") )
    stop( "params must be a numeric matrix of item parameters" )
    
## 3 ## (Make sure that the parameters match the model)
  if( mod == "brm" ){               # for the binary response model:
  
    if( dim(params)[2] != 3 )
      stop( "for a binary response model, there must be three parameters" )
      
    if( any( params[ , 3] >= 1 ) | any( params[ , 3] < 0 ) )
      stop( "the third parameter (guessing) must be between 0 and 1" )
      
    if( any( params[ , 1] < 0 ) )
      warning( "there is at least one negative first parameter (discrimination)" )
      
  } # END if STATEMENTS
  
  if( mod == "grm" ){               # for the graded response model:
  
    if( any( params[ , 1] < 0 ) )
      warning( "there is at least one negative first parameter (discrimination)" )
      
  } # END if STATEMENTS
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Sorting/Naming Parameters #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# If we are using the binary response model, the params need to be named:
  if( mod == "brm" ){
  
    colnames(params) <- c("a", "b", "c")
    
  } # END if STATEMENT

# If we are using the graded response model, the boundaries need to be ordered:
  if( mod == "grm" & inherits(params[ , -1], "matrix") ){
  	
    params           <- cbind(params[ , 1], t(apply(params[ , 2:dim(params)[2] ], MARGIN = 1, FUN = sort)))
    colnames(params) <- c("a", paste("b", 1:(dim(params)[2] - 1), sep = ""))
  
  } # END if STATEMENT
    
    
#~~~~~~~~~~~~~~~~~~~~~~#
# Simulating Responses #
#~~~~~~~~~~~~~~~~~~~~~~#

# Find the response probabilities ("brm" or "grm" or others)
  p <- get( paste0("p.", mod) )(theta, params)
  
# Simulate responses (depending on the model):
  if( mod == "brm" ){
  	sim <- ( p >= runif( length(p) ) ) + 0
  } else{
    sim <- .Call("simPoly", p, ncol(params))
  } # END ifelse STATEMENT
  
# Add the item number as the first column of params:
  params              <- cbind(1:nrow(params), params)
  colnames(params)[1] <- "item"               
    
# Set the class according to the model:
#  a) We want the prime/first class to be "brm", "grm", or whatever model.
#  b) We want the second class to be "matrix" so that it will still ACT like a matrix
#     most of the time except when the model class co-opts.
  class(sim)    <- c(mod, "matrix")
  class(params) <- c(mod, "matrix")

  ret <- list(resp = rbind(sim), params = params, theta = theta)
  
# The above two lines make the code somewhat easy to generalize the model.

  if(mod == "brm" & length(theta) == 1){
    cat("\nBinary response model simulation:\n   ",
         length(theta), " simulee, ",  nrow(params), " items\n\n")
  } else if(mod == "brm"){
     cat("\nBinary response model simulation:\n   ",
         length(theta), " simulees, ", nrow(params), " items\n\n")
  } else if(mod == "grm" & length(theta) == 1){
    cat("\nGraded response model simulation:\n   ",
         length(theta), " simulee, ",  nrow(params), " items, ", ncol(params) - 1, " levels per item\n\n")
  } else if(mod == "grm"){
    cat("\nGraded response model simulation:\n   ",
         length(theta), " simulees, ", nrow(params), " items, ", ncol(params) - 1, " levels per item\n\n")
  } # END ifelse STATEMENTS
                     
  return(ret)
  
} # END simIrt FUNCTION
