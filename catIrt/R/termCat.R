termCat <-
function(params, resp, mod,
         it_flags,
         cat_par, cat_resp, cat_theta,
         cat_info, cat_sem,
         catStart,
         catMiddle,
         catTerm   = list( term  = c("fixed", "precision", "info", "class"),
                           score = c("MLE", "WLE", "BME", "EAP"),
                           n.min = 5, n.max = 50,
                           p.term = list(method   = c("threshold", "change"),
                                         crit     = .25),
                           i.term = list(method   = c("threshold", "change"), 
                                         crit     = 2),
                           c.term = list(method   = c("SPRT", "GLR", "CI"),
                                         bounds   = c(-1, 1),
                                         categ    = c(0, 1, 2),
                                         delta    = .1,
                                         alpha    = .05, beta = .05,
                                         conf.lev = .95)),
         ... )
{
                                     
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# Arguments in catTerm:                                        #
#  - term: how the test should terminate fixed -> var -> class #
#    -- more than one can be selected                          #
#  - n.min, n.max: minimum, maximum number of items for CAT    #
#  - v.term: the variable stopping rule (SEM under ...)        #
#  - c.term: a list of the classification stopping rules       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# NOTE: cat_par need JUST to be the item parameters.

  j <- nrow(cat_par)

#####
# 1 # (CLASSIFICATION DECISION)
#####

# If we need to classify somebody, check classification:
  if( any(catTerm$term == "class") ){
  	
    classFun <- paste("term", catTerm$c.term$method, sep = "")
    categ    <- get(classFun)(cat_par = cat_par, cat_resp = cat_resp,
                              cat_theta = cat_theta, cat_sem = cat_sem,
                              params = params, resp = resp,
                              catMiddle = catMiddle,
                              catTerm = catTerm, ... )
                              
  } else{
  
    categ <- NA
    
  } # END ifelse STATEMENTS


#####
# 2 # (TERMINATION DECISION)
#####

  dec <- stpRule(method    = catTerm$term,       # method used to terminate
                
                 it_min    = catTerm$n.min,      # minimum number of items
                 it_left   = sum(!it_flags),     # number of items left in the bank
                 
                 it_crit   = catTerm$n.max,      # maximum number of items to give
                 se_crit   = catTerm$p.term,     # precision termination criterion
                 info_crit = catTerm$i.term,     # info termination criterion
                 
                 it_obs    = j,                  # number of items given to this point
                 se_obs    = cat_sem,            # standard errors of measurement
                 info_obs  = cat_info,           # observed information
                 
                 categ_est = categ)              # estimated classification to this point    

#####
# 3 # (CLASSIFICATION UNDER TOTAL/FIXED/VARIABLE)
#####

# If they are supposed to be classified, but term was by some other means:
# --> figure out the location of the MLE -- that is the classification.

  if( any(catTerm$term == "class") & ( !is.na(dec$term) & (dec$term != "class") ) ){

    categ <- catTerm$c.term$categ[sum(cat_theta > catTerm$c.term$bounds) + 1]
        
  } # END if STATEMENT
  

  return( list(cat_categ = categ, cat_dec = dec) )
  
} # END termCat FUNCTION
