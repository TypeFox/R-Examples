termGLR <-
function( cat_par, cat_resp, cat_theta,
          catMiddle,
          catTerm, ... )
{

# NOTE 1: cat_par need JUST to be the item parameters.
# NOTE 2: --> "i" is the person index, "j" is the item index, "k" is the bound index.
#         --> "i" was specified in the first (catSim) for loop.

  j      <- length(!is.na(cat_resp))
  c.term <- catTerm$c.term
  likRat <- NULL
  
# Making sure that the category names are what we want:
  if( is.null(c.term$categ) | (length(c.term$categ) != length(c.term$bounds) + 1) )
    c.term$categ <- 1:(length(c.term$bounds) + 1)

# Attach all of the arguments to make them easier to call:

  bounds    <- c.term$bounds # actual bounds
  categ.all <- c.term$categ  # actual categories
  
  delta <- c.term$delta    # half-width of indif
  alpha <- c.term$alpha    # error rate 1
  beta  <- c.term$beta     # error rate 2



#~~~~~~~~~~~~~~~~~~#
# Setting up Stuff #
#~~~~~~~~~~~~~~~~~~#

  c.lower <- log( beta/(1 - alpha) ) # lower critical
  c.upper <- log( (1 - beta)/alpha ) # upper critical

#~~~~~~~~~~~~~#
# Calculating #
#~~~~~~~~~~~~~#

# Hopefully, the class of cat_resp will be "brm" or "grm" at this point.

# The likVals are the likelihood ratio values at points surrounding the cutpoint.
# -> The GLR compares
# ---> a) the likrat at the [thet - delt, thet + delt] with
# ---> b) the likrat at other points.
  theta   <- seq(min(catMiddle$range), max(catMiddle$range), by = .01)
  likVals <- logLik(u = cat_resp, params = cat_par, theta = theta)

# Find the highest point on the likelihood ratio function:
  for( k in seq_along(bounds) ){
  
# Find the likelihood based on the indifference region:
    logLik.u <- max( likVals[theta > bounds[k] + delta] )
    logLik.l <- max( likVals[theta < bounds[k] - delta] )
  
# Standard likelihood ratio ... as always!
    likRat[k] <- logLik.u - logLik.l
      
  } # END for LOOP
  
  likRat <- c(c.upper + .000001, likRat, c.lower - .000001)
  
#~~~~~~~~~#
# Testing #
#~~~~~~~~~#

# Types of decisions: 1, 2, 3, ... , bounds + 1
# No decision: clear = 0
# Decision:    clear = 1

# This is exactly the same as the SPRT:
  for( k in 1:(length(bounds) + 1) ){
  
    if( (likRat[k] >= c.upper) & (likRat[k + 1] <= c.lower) ){
    
      categ <- categ.all[k]

# Make sure to leave, so it doesn't repeat the test.
      break;
    
    } else {
    
      categ <- NA
      
    } # END if STATEMENTS
  
  } # END for LOOP
  
  return( categ )
  
} # END GLR FUNCTION

# Additional?  More complicated bounds based on Bartroff et al. (2008)?
