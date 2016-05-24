middleCat <-
function( params, resp, mod,
          it_flags,
          cat_par, cat_it, cat_resp, cat_theta,
          cat_info, cat_sem,
          catStart,
          catMiddle = list( select = c("UW-FI", "LW-FI", "PW-FI",
                                       "FP-KL", "VP-KL", "FI-KL", "VI-KL",
                                       "random"),
                            at = c("theta", "bounds"),
                            it.range = NULL, n.select = 5,
                            delta = .1,
                            score = c("MLE", "WLE", "BME", "EAP"),
                            range = c(-6, 6),
                            expose = c("none", "SH") ),
          catTerm,
          ddist = dnorm, ... )
{
	
# Make sure that R CMD check doesn't NOTE for binding sake. 
  cat_par.i   <- NULL; rm(cat_par.i)
  cat_it.i    <- NULL; rm(cat_it.i)
  cat_resp.i  <- NULL; rm(cat_resp.i)
  cat_theta.i <- NULL; rm(cat_theta.i)
  cat_info.i  <- NULL; rm(cat_info.i)
  cat_sem.i   <- NULL; rm(cat_sem.i)
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 
# Arguments in catMiddle:                                 #
#  - select: how items will be selected                   #
#  - at: where the items should be selected               #
#  - n.select: number of items to select b/w              #
#  - score: how to score theta after each item            #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  S        <- NULL               # a vector to store selection rates
  j        <- nrow(cat_par) + 1  # the number of items administered (plus 1)
    
  stp <- 0                       # for the S-H iterature item selection

#######################
## I. SELECT AN ITEM ##
#######################

  while(!stp){
    
      it_select <- itChoose( left_par = params[!it_flags, -ncol(params)], mod = mod,
                             numb = catMiddle$n.select, n.select = catMiddle$n.select,
                             cat_par = cat_par[ , -ncol(params)],
                             cat_resp = cat_resp,
                             cat_theta = cat_theta,
                             select = catMiddle$select, at = catMiddle$at,
                             bounds = catTerm$c.term$bounds, delta = catMiddle$delta,
                             range  = catMiddle$range, it.range = catMiddle$it.range,
                             ddist = ddist, ... )$params[ , 1]
                             
# Pick the particular item (using a trick in case we only have one left):
      cat_it.i[j] <<- sample(c(it_select, it_select), size = 1)
    
# Mark the item (getting rid of it), and save the location of the item:
      it_flags[ pl <- which(params[ , 1] == cat_it.i[j]) ] <<- 1
                                              
# IF S-H ITEM EXPOSURE:
#  - a) find the item exposure prob (k),
#  - b) add the item to the total list of selected items,
#  - c) pick a uniform number (u) between 0 and 1
#  - d) If u < k, administer the item, otherwise repeat!
      if(catMiddle$expos == "SH"){
    
        k <- params[pl, ncol(params)]
        S <- c(S, cat_it.i[j])

# If u < k...    	
      if(runif(1) < k){

# ... administer the item and save the parameters:
    	  cat_resp.i[j]  <<- resp[pl]
    	  cat_par.i[j, ] <<- params[pl, ]
    	  stp            <- 1
    	  
    	} # END if STATEMENT
    	
      } else{
      	
# IF NO S-H ITEM EXPOSURE, JUST SAVE THE ITEM:
        S <- c(S, cat_it.i[j])
        cat_resp.i[j]  <<- resp[pl]
        cat_par.i[j, ] <<- params[pl, ]
        stp            <- 1
        
      } # END ifelse STATEMENTS
      
    } # END while STATEMENT
    
#####################
## II. SCORE THETA ##
#####################

## NON-MIXED RESPONSE PATTERN AND MLE ##
  if( { catMiddle$score == "MLE" & 
  	    ( { all(cat_resp.i[1:j] == min(get("resp", envir = environment(middleCat)))) |
  	    	all(cat_resp.i[1:j] == max(get("resp", envir = environment(middleCat)))) } ) } ){
  	    
  	    	
# Repeat the same scoring as in startCat:
    if(catStart$score == "WLE" | catStart$score == "BME" | catStart$score == "EAP"){
    	
# If we are using WLE, EAP, or BME as our method of scoring:
# --> a) Build a character string for the WLE/EAP/BME estimation function,
      scoreFun <- paste(tolower(catStart$score), "Est", sep = "")
      
# --> b) Call the estimation function on stuff to this point,
      x <- get(scoreFun)( resp = cat_resp.i[1:j],
                          params = cat_par.i[1:j, -c(1, ncol(cat_par.i))],
                          range = catMiddle$range, mod = mod, ddist = ddist, ... )
                        
# --> c) Pull out important information.
      cat_theta.i[j + 1] <<- x$theta
      cat_info.i[j]      <<- x$info
      cat_sem.i[j]       <<- x$sem
      
    } else{
    
# If we are doing one of the fixed/random CAT procedures:
# --> a) Figure out which procedure we will use (random/step/fixed),
      cat_theta.i[j + 1] <<- switch(catStart$score,
                                    random = runif(n = 1,
                                                   min = min(catMiddle$range)[1],
                                                   max = max(catMiddle$range)[1]),
                                          
                                    step   = { function(){
                                                if( cat_resp.i[j] == min(get("resp", envir = environment(middleCat))) ){
        	                                      max(cat_theta - catStart$step.size, min(catMiddle$range)[1])
                                                } else if( cat_resp.i[j] == max(get("resp", envir = environment(middleCat))) ){
                                                  min(cat_theta + catStart$step.size, max(catMiddle$range)[1])
                                                } else{
                                                  cat_theta
                                                } # END if STATEMENTS
                                              } }( ),
                                
                                    fixed  = catStart$init.theta, catStart$init.theta)
                                
    } # END ifelse STATEMENTS 
    
    
# After we have given an item and scored it:
# --> BREAK THE FUNCTION AND RETURN STUFF  
  return(list(S = S, j = j))

  } else{
  
# --> a) Build a character string for the WLE/EAP/BME estimation function,
    scoreFun <- paste(tolower(catMiddle$score), "Est", sep = "")
      
# --> b) Call the estimation function on stuff to this point,
    x <- get(scoreFun)( resp = cat_resp.i[1:j],
                        params = cat_par.i[1:j, -c(1, ncol(cat_par.i))],
                        range = catMiddle$range, mod = mod, ddist = ddist, ... )
                        
# --> c) Pull out important information.
    cat_theta.i[j + 1] <<- x$theta
    cat_info.i[j]      <<- x$info
    cat_sem.i[j]       <<- x$sem
      
# After we have given an item and scored it:
# --> BREAK THE FUNCTION AND RETURN STUFF  
    return( list(j = j, S = S) )
      
  } # END ifelse STATEMENTS
  
} # END FUNCTION
