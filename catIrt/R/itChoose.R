itChoose <- function( left_par, mod = c("brm", "grm"),
                      numb = 1, n.select = 1,
                      cat_par = NULL, cat_resp = NULL, cat_theta = NULL,
                      select = c("UW-FI", "LW-FI", "PW-FI",
                                 "FP-KL", "VP-KL", "FI-KL", "VI-KL",
                                 "random"),
                      at = c("theta", "bounds"),
                      range = c(-6, 6), it.range = NULL,
                      delta  = NULL, bounds = NULL,
                      ddist  = dnorm, quad = 33, ... ){                  	

# Note: "at" will be either "theta" or "bounds" depending if the user wants
#        to pick maximum information at "theta" or at one of the classification
#        bounds.
# Note 2: left_par and cat_par need to have and item column.  What to do?
    
# If cat_par and cat_resp don't exist, set them to something.
  if( is.null(cat_par) )
    cat_par  <- c(NA, NA, NA)
  if( is.null(cat_resp) )
    cat_resp <- c(NA, NA, NA)

# Then turn left_par/cat_par into a matrix:
  left_par <- rbind(left_par)
  cat_par  <- rbind(cat_par)
 

#~~~~~~~~~~~~~~~~~#
# Argument Checks #
#~~~~~~~~~~~~~~~~~#

## 1 ## (Make sure that params are ALL numeric)
  if( mode(left_par) != "numeric" )
    stop("left_par needs to be numeric")
 
## 2 ## (Make sure that the dimensions of cat_par and cat_resp match)
  if( all( !is.na(cat_resp) )  & ( length(cat_resp) != dim(cat_par)[1] ) )  
    stop("number of cat_par does not match the length of cat_resp")
    
## 3 ## (If the length of select is strange, set it to FI)
  if( length(select) != 1 )
    select <- "UW-FI"
    
## 4 ## (If the length of mod is strange, pick one)
  if( inherits(left_par, "brm") ){
    mod <- "brm"
  } else if( inherits(left_par, "grm")){
  	mod <- "grm"
  } else if( ( length(mod) != 1 ) | !all( mod %in% c("brm", "grm") ) ){
    stop("mod must be either 'brm' or 'grm'")
  } # END if STATEMENTS
  
## 5 ## (Make sure that delta is OK)
  if( ( length(delta) != 1 ) & ( select %in% c("FP-KL", "VP-KL", "FI-KL", "VI-KL") ) )
    delta <- 1.96
    
#~~~~~~~~~~~~~~~~~#
# Selecting Items #
#~~~~~~~~~~~~~~~~~#

#####
# 1 # (RESTRICTING PARAMETERS)
#####

# If we're using the "brm", b.limits is legit, and we have items within the limits:
#   --> Set left_par to the items within the limits!
  if( mod == "brm" & ( length(it.range) == 2 ) ){
  	
    good_par <- ( left_par[ , 3] >= min(it.range) ) & ( left_par[ , 3] <= max(it.range) )
  
    if( sum(good_par) > 0 )
      left_par <- left_par[good_par, , drop = FALSE]
    
  } # END if STATEMENT
  
#####
# 2 # (WHERE TO SELECT)
#####

  class(left_par) <- class(cat_par) <- c(mod, "matrix")
  
# Next, we need to choose the bound at which we're selecting the item:
  if( (length(at) != 1) | all(at == "theta") | is.null(bounds) ){
  
    theta <- cat_theta
    
  } else if( at == "bounds" ){
  	
# Figuring out the closest bound to current theta, and sampling lest there are 2:
    theta  <- bounds[ which.min( abs( cat_theta - bounds ) ) ]
  
  }
  
# Note: We might want more sophsticated methods of choosing proximate theta?

#####
# 3 # (FIGURING OUT INFO)
#####

## FOR FISHER INFORMATIONS ##
#  -- we need expected info for remaining params at current theta:
#  -- we need to figure out which boundary we should be using
#  -- we need to figure out Fisher Info at that boundary
  if( { select == "UW-FI" |
  	    ( ( select %in% c("LW-FI", "PW-FI") ) & ( any( is.na(cat_resp) ) | ( length(cat_resp) == 0 ) ) ) } ){
  
    info <- FI(params = left_par[ , -1], theta = theta, type = "expected")$item
    
  } else if( select == "LW-FI" ){
  
    info <- WFI( left_par = left_par[ , -1], mod = mod,
                 cat_par = cat_par[ , -1], cat_resp = cat_resp,
                 range = range,
                 type = "likelihood", quad = quad )$item
    
  } else if( select == "PW-FI" ){
  	
    info <- WFI( left_par = left_par[ , -1], mod = mod,
                 cat_par = cat_par[ , -1], cat_resp = cat_resp,
                 range = range,
                 type = "posterior", quad = quad,
                 ddist = ddist, ... )$item
    
  } # END ifelse STATEMENTS
  
# Likelihood Weighted and Posterior Weighted are calculated differently.


## FOR KL INFORMATION ##
#  -- we need to figure out which boundary we should be using,
#  -- we need to find KL Info at that boundary point using delt,
#  -- if variable KL, we correct using the sqrt of the number of items,
#  -- if integrated KL, we integrate from theta - delta to theta + delta.
  if( select == "FP-KL" ){
      
    info <-  KL(params = left_par[ , -1],
                theta  = theta,
                delta  = delta / 1)$item
    
  } else if( select == "VP-KL" ){
  	
  	info <-  KL(params = left_par[ , -1],
  	            theta  = theta,
  	            delta  = delta / sqrt( max(nrow(cat_par), 1) ) )$item
  	
  }	else if( select == "FI-KL" ){
  
    info <- IKL(params = left_par[ , -1],
                theta  = theta,
                delta  = delta / 1,
                quad   = quad)$item
  
  } else if( select == "VI-KL" ){
  	
    info <- IKL(params = left_par[ , -1], 
                theta  = theta,
                delta  = delta / sqrt( max(nrow(cat_par), 1) ),
                quad   = quad)$item
    
  } # END ifelse STATEMENTS
  
  
## FOR random INFORMATION ##
#  -- everything should have the same information
#  -- or at least, probabilistically, everything should have the same information :)
#  -- the runif is to make sure that the we don't just select items at the beginning.
  if( select == "random" ){

    info <- runif( nrow(left_par) )
  
  } # END random if STATEMENT

#####
# 4 # (SORTING ITEMS)
#####

# Making sure that we have items to choose between:
  n_select <- ifelse(test = is.null(n.select), yes = 1, no = n.select)
  
# In case there are fewer items left than those to choose from:
  n_select <- min(n_select, length(info))[1]  # for the items to choose between
  numb     <- min(numb, length(info))[1]      # for the number of items to select

# Next we want to find the "n" items with largest info (making sure to keep everything as a matrix):
  par_select  <- left_par[order(info, decreasing = TRUE), , drop = FALSE][1:max(numb, n_select)[1], , drop = FALSE]
  info_select <- sort(info, decreasing = TRUE)[1:max(numb, n_select)[1]]
  
#####
# 5 # (CHOOSING ITEMS)
#####

# If we are only selecting one item, randomly select that item.
  if( numb == 1 ){

# (doubling the number of items is a trick to make sure it selects at least one):
    it_next   <- sample(c( par_select[ , 1], par_select[ , 1]), size = 1)
    par_next  <- par_select[par_select[ , 1] == it_next, , drop = FALSE]
    info_next <- info_select[par_select[ , 1] == it_next]
    
  } else if( (numb > 1) & (numb < n_select) ){
  
# If we are selecting more than one item, randomly select those items:
    it_next   <- sample( c(par_select[ , 1]), size = numb)
    par_next  <- par_select[par_select[ , 1] %in% it_next, , drop = FALSE]
    info_next <- info_select[par_select[ , 1] %in% it_next]
    
  } else if( numb >= n_select ){
  
# If we are selecting all of the items we have, select all of them:
    par_next  <- par_select
    info_next <- info_select
    
  } # END if STATEMENTS
  
  
  return( list(params = par_next, info = info_next, type = select) )
  
} # END itChoose FUNCTION