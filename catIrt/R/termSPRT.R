termSPRT <-
function( cat_par, cat_resp, cat_theta,
          catTerm, ... ){
          	
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

# The logLikelihood ratio for each bounds
  for( k in seq_along(bounds) ){
  
    likRat[k] <- { logLik(u = cat_resp, params = cat_par, theta = bounds[k] + delta) -
                   logLik(u = cat_resp, params = cat_par, theta = bounds[k] - delta) }
      
  } # END for LOOP
  
  likRat <- c(c.upper + .000001, likRat, c.lower - .000001)
  
#~~~~~~~~~#
# Testing #
#~~~~~~~~~#

# Types of decisions: 1, 2, 3, ... , bounds + 1
# No decision: clear = 0
# Decision:    clear = 1

  for( k in 1:(length(bounds) + 1) ){
  
    if( (likRat[k] >= c.upper) & (likRat[k + 1] <= c.lower) ){
    
      categ <- categ.all[k]
      clear <- 1

      break
    
    } else {
    
      categ <- NA
      clear <- 0
      
    } # END if STATEMENTS
  
  } # END for LOOP
  
  return( categ )
  
} # END SPRT FUNCTION

# Additional?  Change calculation or critical values based on Spray (1993)?
