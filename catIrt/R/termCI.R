termCI <-
function( cat_par, cat_resp,
          cat_theta, cat_sem,
          catTerm, ... ){
          	
# NOTE 1: cat_par need JUST to be the item parameters.
# NOTE 2: --> "i" is the person index, "j" is the item index, "k" is the bound index.
#         --> "i" was specified in the first (catSim) for loop.
                 
  j      <- length(!is.na(cat_resp))
  c.term <- catTerm$c.term

# Making sure that the category names are what we want:
  if( is.null(c.term$categ) | (length(c.term$categ) != length(c.term$bounds) + 1) )
    c.term$categ <- 1:(length(c.term$bounds) + 1)

# Attach all of the arguments to make them easier to call:

  bounds    <- c.term$bounds # actual bounds
  categ.all <- c.term$categ  # actual categories
  
#   1) Calculate estimated theta,
#   2) Calculated SEM based on items taken,
#   3) Calculate Confidence Int: theta-hat +/- 1.96*SEM
#   4) If theta-hat +/- 1.96*SEM is within bounds, take that level.

#~~~~~~~~~~~~~#
# Calculating #
#~~~~~~~~~~~~~#

# First, calculate estimated theta (cat_thet)

# Second, calculate SEM based on items taken:
  std_err <- cat_sem[j]
  
  if( is.na(std_err) )
    std_err <- Inf
  
# Third, calculating Confidence Interval stuff:
  theta <- cat_theta + c(-1, 1)*qnorm( (1 + c.term$conf.lev)/2 ) * std_err
  
    
# Making sure we can check every condition!
  bounds <- c(-Inf, bounds, Inf)

#~~~~~~~~~#
# Testing #
#~~~~~~~~~#

# Types of decisions: 1, 2, 3, ... , bounds + 1
# No decision: clear = 0
# Decision:    clear = 1

  for( k in 1:(length(bounds) - 1) ){
  
    if( theta[1] > bounds[k] & theta[2] < bounds[k + 1] ){
    
      categ <- categ.all[k]

# Make sure to leave, so it doesn't repeat the test.
      break;
    
    } else {
    
      categ <- NA
      
    } # END if STATEMENTS
  
  } # END for LOOP
  
  return( categ )
  
} # END CI FUNCTION

# Additional?  Correction for multiple testing?

