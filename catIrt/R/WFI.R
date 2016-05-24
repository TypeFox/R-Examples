WFI <- function( left_par, mod = c("brm", "grm"),
                 cat_par, cat_resp,
                 range = c(-6, 6),
                 type = c("likelihood", "posterior"),
                 quad = 33,
                 ddist = dnorm, ... )
{
  
# This function finds Weighted Fisher Information:
# -- Weighted on the likelihood function, or
# -- Weighted on the posterior information function.

# Both situations require:
# a) The parameters to this point in the CAT, and
# b) The responses to this point in the CAT.

# Turn params into a matrix:
  left_par <- rbind(left_par)
  cat_par  <- rbind(cat_par)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Information and Likelihood #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Indicate the lower/upper bounds of integration and density (if it doesn't exist):
  if( is.null(range) )
    range <- c(-6, 6)
  if( is.null(ddist) )
    ddist <- dnorm
    
  l <- range[1]; u <- range[2]
  
# Make sure that the lower/upper points of integration is within the legal range:
  tmp.x <- seq(l, u, by = .01)
  tmp.y <- ddist(x = tmp.x, ... )
  
# Usually, if the density doesn't exist, it prints "0", but it might print NA:
  l <- min( tmp.x[signif(tmp.y) != 0 & !is.na(tmp.y)] )
  u <- max( tmp.x[signif(tmp.y) != 0 & !is.na(tmp.y)] )
  
# Set up the bounds of integration and turn cat_resp into a matrix:
  X        <- seq(l, u, length = quad)
  cat_resp <- matrix(cat_resp, nrow = quad, ncol = nrow(cat_par), byrow = TRUE)
  
# The Likelihood function used in the integration:
  LikFun  <- function( ... )
    exp( get(paste("logLik.", mod, sep = "") )( ... ) )
  FishFun <- function( ... )
    get( paste("FI.", mod, sep = "") )( ... )
               
# The Likelihood and FI for all of the items thus far:
  wfi.lik  <- LikFun(u = cat_resp, theta = X, params = cat_par)
  wfi.fish <- FishFun(params = left_par, theta = X, type = "expected")$item  
  
  if( nrow(left_par) == 1 )
    wfi.fish <- cbind(wfi.fish)
               
# Note: Turning cat_resp into a matrix is probably inefficient but transparent.

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Set-Up Integral Function #
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

  if( type == "likelihood" ){
  	
    Y <- wfi.fish * wfi.lik
    
  } else if( type == "posterior" ){
    
    Y <- wfi.fish * wfi.lik * ddist(x = X, ... ) 
   
  } # END ifelse FUNCTIONS

#~~~~~~~~~~~~~~~#
# And Integrate #
#~~~~~~~~~~~~~~~#

  info <- apply(Y, MARGIN = 2, FUN = integrate.q, x = X)
  
  return( list(item = info, type = type) )
  
} # END WFI FUNCTION
