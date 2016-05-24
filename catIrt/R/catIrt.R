# This function will be the main CAT function.

catIrt <- function( params, mod = c("brm", "grm"),
                    resp      = NULL,
                    theta     = NULL,
                    catStart  = list( n.start = 5, init.theta = 0,
                                      select = c("UW-FI", "LW-FI", "PW-FI",
                                                 "FP-KL", "VP-KL", "FI-KL", "VI-KL",
                                                 "random"),
                                      at = c("theta", "bounds"),
                                      it.range = NULL, n.select = 1,
                                      delta = .1,
                                      score = c("fixed", "step", "random",
                                                "WLE", "BME", "EAP"),
                                      range = c(-1, 1),
                                      step.size = 3, leave.after.MLE = FALSE ),
                    catMiddle = list( select = c("UW-FI", "LW-FI", "PW-FI",
                                                 "FP-KL", "VP-KL", "FI-KL", "VI-KL",
                                                 "random"),
                                      at = c("theta", "bounds"),
                                      it.range = NULL, n.select = 1,
                                      delta = .1,
                                      score = c("MLE", "WLE", "BME", "EAP"),
                                      range = c(-6, 6),
                                      expos = c("none", "SH") ),
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
                    ddist     = dnorm,
                    progress  = TRUE, ... )                                  
{

# Make sure that the environments are OK, so that the "<<-" works.
  environment(startCat)  <- environment()
  environment(middleCat) <- environment()
  environment(termCat)   <- environment()


#############################################################################
######################## BEGIN ARGUMENT CHECK SECTION #######################

#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ARGUMENT CHECKS (PART 1) # (INITIAL CHECKING OF: resp/params/theta/mod)
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

## 1 ## (Make sure that the parameter matrix exists)
  if( missing(params) )
    stop( "need to supply a 'params' matrix" )


## 2 ## (Make sure that one of "resp" or "thetas" exists and is of appropriate dimension)
  if( is.null(resp) & is.null(theta) ){
    stop( "need to include 'resp' and/or 'theta'" )
    
  } else if( !is.null(resp) & all(dim(resp)[2] != dim(params)[1]) ){
    stop( "the number of 'resp' columns must equal the number of 'params' rows" )
    
  } else if( ( !is.null(resp) & !is.null(theta) ) & all(dim(resp)[1] != length(theta)) ){
  	warning( "the number of 'resp' rows does not equal the length of 'theta'" )

  } # END ifelse STATEMENTS

  
## 3 ## (Make sure that the model is appropriately declared)
  mod.opt <- c("brm", "grm")

  if(missing(mod))
    mod  <- NULL
    
  if( (length(mod) != 1) | !any(mod %in% mod.opt) ){
  	
# --> If class of resp inherits any of the appropriate models, set mod to that model OR
    if( any(class(resp) %in% mod.opt) ){
      mod <- class(resp)[class(resp) %in% mod.opt][1]  
        	  
# --> Ask the user for the model (repeating until you get it).
    } else if( interactive() ){       			
      while( !( length(mod) == 1 & all(mod %in% mod.opt) ) )
        mod <- readline( paste("Select from ONLY ONE of the following IRT models - ",
    	                       paste(mod.opt, collapse = ", "),
    	                       ": ", sep = "")
    	               )    	                 	                    
    } else{
        stop( paste("'mod' must be ONLY ONE of ", paste(mod.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
  } # END if STATEMENT


## 4 ## (If any of catStart/catMiddle/catTerm doesn't exist, then initialize it)
  if(missing(catStart))
    catStart  <- NULL
  if(missing(catMiddle))
    catMiddle <- NULL
  if(missing(catTerm))
    catTerm   <- NULL


## 5 ## (Make sure the parameters are OK for exposure control)

# --> If the user specifies multiple exposures or not an appropriate one, exposure should be none.
  if( (length(catMiddle$expos) != 1) | !any(catMiddle$expos %in% c("none", "SH")) )
    catMiddle$expos <- "none"	
  	   
# --> If exposure is none, add an P(S | A) column of 1s.
  if( catMiddle$expos == "none" )
    params <- cbind(params, 1)
  
# --> If the user exists, warn the user about setting up the param matrix appropriately.  
  if( interactive() & progress )
    cat("\nIff 'catMiddle$expos' equals 'SH', the last 'params' column must contain P(Admin | Select)'s.\n\n")    


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# BUILDING THE PARAMS/RESPONSE MATRIX #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Making sure the parameters are of a particular class:
  if( inherits(params, mod.opt) ){
    if( any( round(params[ , 1]) != 1:nrow(params) ) ){
      params <- cbind(1:nrow(params), params)  # first column is item number
    } # END if STATEMENT
  } else{
    params <- cbind(1:nrow(params), params)    # first column is item number
  } # END ifelse STATEMENT

# Building the response matrix and indicating its class:
  if( is.null(resp) ){
    resp <- simIrt(theta = theta, params = params[ , -c(1, ncol(params))], mod = mod)$resp
    
  } else{
  	
# And then make sure that 'resp' is a matrix:
  	resp        <- rbind(resp)
  	class(resp) <- c(mod, "matrix")
  	
  } # END ifelse STATEMENT
  
# And then make sure that 'resp' is a matrix:
  if( is.null( dim(resp) ) ){                # if it's a vector ... -->
    resp        <- matrix(resp, nrow = 1)    # ... --> turn it into a multi-column matrix,
    class(resp) <- c(mod, "matrix")          # ... --> and indicate its class
  }  # END if STATEMENT
 
   
#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ARGUMENT CHECKS (PART 2) # (MODEL-SPECIFIC CHECKS)
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

## FOR THE BINARY RESPONSE MODEL ##
  if( mod == "brm" ){

## 1 ## (Make sure resp has all 0's and 1's and warn if it's ALL 0's or 1's)
    if( any(resp != 0  & resp != 1) )
      stop( "'resp' must only contain 0's and 1's" )
  
    if( all(resp == 0) | all(resp == 1) )
      warning( "'resp' contains only 0's or only 1's" )

## 2 ## (Make sure that the dimensions of params are correct)
    if( dim(params)[2] != 5 )
      stop( "'params' must have three or four columns: item, a, b, and c" )

## 3 ## (Naming and classing the parameters)
    colnames(params) <- c("item", "a", "b", "c", "SH")
    class(params)    <- c(mod, "matrix") 
      
  } # END if (brm) STATEMENT


## FOR THE GRADED RESPONSE MODEL ##
  if( mod == "grm" ){
    
## 1 ## (Make sure that the responses are all positive)
    if( any(resp < 0) )
      stop( "'resp' must contain positive integers" )
      
## 2 ## (Make sure that the responses are not too big)
    if( any(resp > dim(params)[2]) )
      stop( "every entry of 'resp' must be less than the number of columns of 'params'" )
      
## 3 ## (Make sure the parameters are in the appropriate order)
    params <- cbind( params[, 1:2],
                     t(apply(params[, 3:(dim(params)[2] - 1)], MARGIN = 1, FUN = sort)),
                     params[, dim(params)[2]] )
                    
## 4 ## (Naming and classing the parameters)
    colnames(params) <- c("item", "a", paste("b", 1:(dim(params)[2] - 3), sep = ""), "SH")
    class(params)    <- c(mod, "matrix") 
      
  } # END if (grm) STATEMENT


#~~~~~~~~~~~~~~~~~~~~~~~~~~#
# ARGUMENT CHECKS (PART 3) # (CHECKING OF: catStart/catMiddle/catTerm)
#~~~~~~~~~~~~~~~~~~~~~~~~~~#

## I. FOR THE CATSTART LIST ##
	
## LIST OF OPTIONS ##
  sel.opt <- c("UW-FI", "LW-FI", "PW-FI", "FP-KL", "VP-KL", "FI-KL", "VI-KL", "random")
  at.opt  <- c("theta", "bounds")
  sco.opt <- c("fixed", "step", "random", "WLE", "BME", "EAP")

## a) n.start ##     
  if( {length(catStart$n.start) != 1 | !all(catStart$n.start %in% 1:nrow(params)) } ){

# --> Make sure 'n.it' is a positive integer that is not too large of a number.
    if( interactive() ){
      while( !( length(catStart$n.start) == 1 & all(catStart$n.start %in% 1:nrow(params)) ) ){
        catStart$n.start <- readline("Select the number of items for starting the CAT: ")
        catStart$n.start <- suppressWarnings(as.numeric(catStart$n.start))         
      } # END while LOOP
    } else{
      catStart$n.start <- 1
    } # END ifelse STATEMENTS
  
  } # END if STATEMENT
  
  
## b) init.theta ## (a number that ends up being as long as the number of examinees)

# --> If init.theta isn't specified, then it should be NA so the functions work.
  if( is.null(catStart$init.theta) )
    catStart$init.theta <- NA
    
# --> If init.theta is an odd length, warn the user.
  if( length(catStart$init.theta) != 1 & (length(catStart$init.theta) != dim(resp)[1]) )
    warning( "initial values are not specified for each simulee" )

# --> Turn init.theta into a numeric variable.    
  catStart$init.theta <- suppressWarnings(as.numeric(catStart$init.theta))
  
# --> And make sure that it is an actual number (and NOT NA).
  if( interactive() ){
  	while( any(is.na(catStart$init.theta)) ){
      catStart$init.theta <- readline("Select the initial/fixed trait estimate used for all simulees: ")
  	  catStart$init.theta <- suppressWarnings(as.numeric(catStart$init.theta))
  	} # END while LOOP 
  } else{
    catStart$init.theta <- rep(0, length.out = dim(resp)) 
  } # END ifelse STATEMENT

# Finally, we want the same number of initial ability values as there are simulees.  
  catStart$init.theta <- rep(catStart$init.theta, length.out = dim(resp)[1])


## c) select ##
  if( (length(catStart$select) != 1) | !any(catStart$select %in% sel.opt) ){
  	
# --> Make sure 'select' matches one of the possible selection mechanisms.
    if( interactive() ){    			
      while( !( length(catStart$select) == 1 & all(catStart$select %in% sel.opt) ) )
        catStart$select <- readline( paste("Select from ONLY ONE of the following starting methods to select items - ",
    	                                   paste(sel.opt, collapse = ", "),
    	                                   ": ", sep = "")
    	                           )    	                 	                    
    } else{
      stop( paste("'catStart$select' must be ONLY ONE of ", paste(sel.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
    
  } # END if STATEMENT


## d) at ## (only if 'select' is not 'random', 'LW-FI' ,or 'PW-FI')
  if( (length(catStart$at) != 1) | !any(catStart$at %in% at.opt) ){
 
# --> If we are randomly selecting or selecting based on weights, who cares about where.	
  	if( catStart$select == "random" | catStart$select == "LW-FI" | catStart$select == "PW-FI" ){
  	  catStart$at <- "theta"

# --> Make sure 'at' matches one of the possible selection-at mechanisms.
  	} else if( interactive() ){       			
      while( !( length(catStart$at) == 1 & all(catStart$at %in% at.opt) ) )
        catStart$at <- readline( paste("Select from ONLY ONE of the following starting locations to select items - ",
    	                               paste(at.opt, collapse = ", "),
    	                               ": ", sep = "")
    	                       )    	                 	                    
    } else{
      stop( paste("'catStart$at' must be ONLY ONE of ", paste(at.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
    
  } # END if STATEMENT
  

## e) it.range ##
# --> If 'it.range' is specified for non-binary response models, warn
  if( mod != "brm" & !is.null(catStart$it.range) )
    warning("'catStart$it.range' can only be specified for 'brm'")

## f) n.select ##

# --> If 'n.select' is not specified, set it to 1.
  if( length(catStart$n.select) != 1 )
    catStart$n.select <- 1
  
  
## g) delta ## (only if using KL selection mechanism)
  if( any( catStart$select %in% c("FP-KL", "VP-KL", "FI-KL", "VI-KL") ) ){
  	
    if( length(catStart$delta) != 1 )
      catStart$delta <- NA
    
    if( any(is.na(catStart$delta)) | any(!is.numeric(catStart$delta)) ){
 
# --> Make sure 'delta' is a positive number.    	
      if(interactive()){
        while( !( { length(catStart$delta) == 1 &
        	            all(!is.na(catStart$delta)) &
        	            all(is.numeric(catStart$delta)) } ) ){
          catStart$delta <- readline("Select a starting half-width constant for use in KL information: ")
          catStart$delta <- suppressWarnings(abs(as.numeric(catStart$delta)))
        } # END while LOOP
      } else{
        catStart$delta <- .1
      } # END ifelse STATEMENTS
       
    } # END if STATEMENT  
    
  } # END if STATEMENT


## h) score ##
  if( (length(catStart$score) != 1) | !any(catStart$score %in% sco.opt) ){ 

# --> Make sure 'score' matches one of the possible scoring methods.   
    if( interactive() ){    			
      while( !( length(catStart$score) == 1 & all(catStart$score %in% sco.opt) ) )
        catStart$score <- readline( paste("Select from ONLY ONE of the following starting methods to score the latent trait - ",
    	                                  paste(sco.opt, collapse = ", "),
    	                                  ": ", sep = "")
    	                          )     	                 	                    
    } else{
      stop( paste("'catStart$score' must be ONLY ONE of ", paste(sco.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
    
  } # END if STATEMENT


## i) range ## (Make sure that the MLE/EAP/BME has an integer to maximize)

# --> If 'int' is not correctly specified, set it to a default.
  if( length(catStart$range) != 2 )
    catStart$range <- c(-6, 6)
    
# --> And making sure that the interval is in sorted order.
  catStart$range <- sort(catStart$range)
  
    
## j) step.size ## (only if 'score' is 'step')
  if( { catStart$score == "step" &
  	    ( length(catStart$step.size) != 1 | !is.numeric(catStart$step.size) ) } ){

# --> Make sure 'step.size' is a positive number. 	    	 
    if( interactive() ){
      catStart$step.size <- NA

# --> If we don't have a number yet, read in a value and try to make it a positive number.      			
      while( !( { length(catStart$step.size) == 1 &
      	          all( !is.na(catStart$step.size) ) & is.numeric(catStart$step.size) } ) ){
      	catStart$step.size <- readline("Select the initial change in latent trait for a correct or incorrect response: ")           	
        catStart$step.size <- suppressWarnings(abs(as.numeric(catStart$step.size)))       
      } # END while STATEMENT    	                 	                    
    } else{
      catStart$step.size <- 1
    } # END ifelse STATEMENTS
  } # END if STATEMENT
	
  	  
## k) leave.after.MLE ##

# --> If 'leave.after.MLE' is not specified, set it to FALSE.
  if( !is.logical(catStart$leave.after.MLE) )
    catStart$leave.after.MLE <- FALSE


## II. FOR THE CATMIDDLE LIST ##

## LIST OF OPTIONS ##
  sel.opt <- c("UW-FI", "LW-FI", "PW-FI", "FP-KL", "VP-KL", "FI-KL", "VI-KL", "random")
  at.opt   <- c("theta", "bounds")
  sco.opt  <- c("MLE", "WLE", "BME", "EAP")


## a) select ##
  if( (length(catMiddle$select) != 1) | !any(catMiddle$select %in% sel.opt) ){
 
# --> Make sure 'select' matches one of the possible selection mechanisms.
    if( interactive() ){    
      while( !( length(catMiddle$select) == 1 & all(catMiddle$select %in% sel.opt) ) )
        catMiddle$select <- readline( paste("Select from ONLY ONE of the following middle methods to select items - ",
    	                                    paste(sel.opt, collapse = ", "),
    	                                    ": ", sep = "")
    	                            )    	                 	                    
    } else{
      stop( paste("'catMiddle$select' must be ONLY ONE of ", paste(sel.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
  } # END if STATEMENT


## b) at ## (only if 'select' is not 'random', 'LW-FI' ,or 'PW-FI')
  if( (length(catMiddle$at) != 1) | !any(catMiddle$at %in% at.opt) ){
 
# --> If we are randomly selecting or selecting based on weights, who cares about where.	
  	if( catMiddle$select == "random" | catMiddle$select == "LW-FI" | catMiddle$select == "PW-FI" ){
  	  catMiddle$at <- "theta"

# --> Make sure 'at' matches one of the possible selection-at mechanisms.
  	} else if( interactive() ){       			
      while( !( length(catMiddle$at) == 1 & all(catMiddle$at %in% at.opt) ) )
        catMiddle$at <- readline( paste("Select from ONLY ONE of the following middle locations to select items - ",
    	                               paste(at.opt, collapse = ", "),
    	                               ": ", sep = "")
    	                       )    	                 	                    
    } else{
      stop( paste("'catStart$at' must be ONLY ONE of ", paste(at.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
    
  } # END if STATEMENT
  
## c) it.range ##
# --> If 'it.range' is specified for non-binary response models, warn
  if( mod != "brm" & !is.null(catMiddle$it.range) )
    warning("'catMiddle$it.range' can only be specified for 'brm'")
    
## d) n.select ##

# --> If 'n.select' is not specified, set it to 1.
  if( length(catMiddle$n.select) != 1)
    catMiddle$n.select <- 1
  
## e) delta ## (only if using KL selection mechanism)
  if( any( catMiddle$select %in% c("FP-KL", "VP-KL", "FI-KL", "VI-KL") ) ){
  	
    if( length(catMiddle$delta) != 1 )
      catMiddle$delta <- NA
    
    if( any(is.na(catMiddle$delta)) | any(!is.numeric(catMiddle$delta)) ){
 
# --> Make sure 'delta' is a positive number.    	
      if(interactive()){
        while( !( { length(catMiddle$delta) == 1 &
        	            all(!is.na(catMiddle$delta)) &
        	            all(is.numeric(catMiddle$delta)) } ) ){
          catMiddle$delta <- readline("Select a middle half-width constant for use in KL information: ")
          catMiddle$delta <- suppressWarnings(abs(as.numeric(catMiddle$delta)))
        } # END while LOOP
      } else{
        catMiddle$delta <- .1
      } # END ifelse STATEMENTS
       
    } # END if STATEMENT  
    
  } # END if STATEMENT


## f) score ##
  if( (length(catMiddle$score) != 1) | !any(catMiddle$score %in% sco.opt) ){
  	
# --> Make sure 'score' matches one of the possible scoring methods.       
    if( interactive() ){        			
      while( !( length(catMiddle$score) == 1 & all(catMiddle$score %in% sco.opt) ) )
        catMiddle$score <- readline( paste("Select from ONLY ONE of the following middle methods to score the latent trait - ",
    	                                   paste(sco.opt, collapse = ", "),
    	                                   ": ", sep = "")
    	                           )    	                 	                    
    } else{
      stop( paste("'catMiddle$score' must be ONLY ONE of ", paste(sco.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
  } # END if STATEMENT


## g) range ## (Make sure that the MLE/EAP/BME has an integer to maximize)

# --> If 'int' is not correctly specified, set it to a default.
  if( length(catMiddle$range) != 2 )
    catMiddle$range <- c(-6, 6)
    
# --> And making sure that the interval is in sorted order.
  catMiddle$range <- sort(catMiddle$range)
  
  
## III. FOR THE CATTERM LIST ##

## LIST OF OPTIONS ##
  term.opt   <- c("fixed", "precision", "info", "class")
  sco.opt    <- c("MLE", "WLE", "BME", "EAP")
  p.meth.opt <- c("threshold", "change")
  i.meth.opt <- c("threshold", "change")
  c.meth.opt <- c("SPRT", "GLR", "CI")
  

  
## a) term ##
  if( (length(catTerm$term) < 1) | !all(catTerm$term %in% term.opt) ){

# --> Make sure 'term' matches at least one of the possible termination criteria.   	  
    if( interactive() ){        			
      while( !( length(catTerm$term) >= 1 & all(catTerm$term %in% term.opt) ) ){
        cat("Select from one or more of the following stopping rules (press 'return' after each) - ",
            paste(term.opt, collapse = ", "),
            ": ", sep = "")
        catTerm$term <- scan(what = character(), n = length(term.opt), quiet = TRUE)
      } # END while STATEMENT                  	                 	                    
    } else{
      stop( paste("'catTerm$term' must be ONLY ", paste(term.opt, collapse = " and/or "), sep = "" ) )
    } # END ifelse STATEMENTS
    
  } # END if STATEMENT


## b) score ##

# --> If score at the end of the cat isn't specified, set it to the typical scoring.
  if( (length(catTerm$score) < 1) | !all(catTerm$score %in% sco.opt))
    catTerm$score <- catMiddle$score


## c) n.min ## (must be an integer greater than or equal to 0 and less than the number of params)
  if( { is.null(catTerm$n.min) |
  	    !all(catTerm$n.min %in%  1:nrow(params)) } ){

# --> Make sure 'n.it' is a positive integer that is not too large of a number. 	   
    if( interactive() ){    			
      while( !( length(catTerm$n.min) == 1 & all(catTerm$n.min %in%  1:nrow(params)) ) ){
        catTerm$n.min <- readline("Select a minimum number of items: ")
        catTerm$n.min <- suppressWarnings(as.numeric(catTerm$n.min))  
      } # END while LOOP
    } else{
      catTerm$n.min <- 1
    } # END ifelse STATEMENTS
  
  } # END if STATEMENT

 
## d) n.max ## (must be an integer greater than n.min AND n.it or "all")
  if( { is.null(catTerm$n.max) |
  	    !all(catTerm$n.max %in%  max(catTerm$n.min, catStart$n.it):nrow(params)) } ){

    if(catTerm$n.max == "all" | catTerm$n.max == "'all'" | catTerm$n.max == "\"all\"")
      catTerm$n.max <- nrow(params)
      
# --> Make sure 'n.it' is a positive integer that is not too small of a number. 	    	  
    if( interactive() ){			
      while( !( { (length(catTerm$n.max) == 1) &
      	          all(catTerm$n.max %in%  max(catTerm$n.min, catStart$n.it):nrow(params)) } ) ){      	
        catTerm$n.max <- readline("Select a maximum number of items or 'all': ")

# --> If "all" is selected, break the while loop so the comparison isn't made.
        if(catTerm$n.max == "all" | catTerm$n.max == "'all'" | catTerm$n.max == "\"all\""){
          catTerm$n.max <- nrow(params)
        } else{
          catTerm$n.max <- suppressWarnings(as.numeric(catTerm$n.max))
        } # END ifelse STATEMENT
     
      } # END while LOOP
    	                 	                    
    } else{
      catTerm$n.max <- nrow(params)
    } # END ifelse STATEMENTS

  } # END if STATEMENT
  
# And if n.max is less than the total number of params, one of the termination is implicitly fixed.
  if( catTerm$n.max < nrow(params) )
    catTerm$term <- c(catTerm$term, "fixed")
    
  
## e) p.term ##

## e1) p.term --> method ## (only for the precision termination criterion)
  if( { any(catTerm$term == "precision") &
  	    ( (length(catTerm$p.term$method) != 1) | !any(catTerm$p.term$method %in% p.meth.opt) ) } ){

# --> Make sure 'method' matches at least one of the possible method stopping rules.
    if( interactive() ){
      while( !( length(catTerm$p.term$method) == 1 & all(catTerm$p.term$method %in% p.meth.opt) ) )
        catTerm$p.term$method <- readline( paste("Select from ONLY ONE of the following precision-based stopping rules - ",
    	                                         paste(p.meth.opt, collapse = ", "),
    	                                         ": ", sep = "")
    	                                 )    	                 	                    
    } else{
      stop( paste("'catTerm$p.term$method' must be ONLY ONE of ", paste(p.meth.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
    
  } # END if STATEMENT

  
## e2) p.term --> crit ##

## threshold ##
  if( any(catTerm$term == "precision" & catTerm$p.term$method %in% c("threshold") ) ){
  
    if( length(catTerm$p.term$crit) != 1 )
      catTerm$p.term$crit <- NA
      
    if( any(is.na(catTerm$p.term$crit)) | any(!is.numeric(catTerm$p.term$crit)) ){
    	
# --> Make sure 'crit' is a positive number.
      if( interactive() ){
        repeat{
          catTerm$p.term$crit <- readline("Select an SEM cut-off for variable termination: ")
          catTerm$p.term$crit <- suppressWarnings(as.numeric(catTerm$p.term$crit))
          
          if( !( all(!is.na(catTerm$p.term$crit)) & all(is.numeric(catTerm$p.term$crit)) ) )
            next;
          
          if( catTerm$p.term$crit >= 0 )
            break;
            
        } # END repeat LOOP
      } else{
        catTerm$p.term$crit <- .25
      } # END ifelse STATEMENTS
    } # END if STATEMENT

  } # END if STATEMENT


## change ##

  if( any(catTerm$term == "precision" & catTerm$p.term$method %in% c("change") ) ){
  
    if( length(catTerm$p.term$crit) != 1 )
      catTerm$p.term$crit <- NA
      
    if( any(is.na(catTerm$p.term$crit)) | any(!is.numeric(catTerm$p.term$crit)) ){
    	
# --> Make sure 'crit' is a positive number.
      if( interactive() ){
        repeat{
          catTerm$p.term$crit <- readline("Select the maximum SEM change for variable termination: ")
          catTerm$p.term$crit <- suppressWarnings(as.numeric(catTerm$p.term$crit))
          
          if( !( all(!is.na(catTerm$p.term$crit)) & all(is.numeric(catTerm$p.term$crit)) ) )
            next;
          
          if( catTerm$p.term$crit >= 0 )
            break;
            
        } # END repeat LOOP
      } else{
        catTerm$p.term$crit <- .25
      } # END ifelse STATEMENTS
    } # END if STATEMENT

  } # END if STATEMENT

  
## f) i.term ## (only for the info termination criterion)

## f1) i.term --> method ## (only for the info termination criterion)
  if( { any(catTerm$term == "info") &
  	    ( (length(catTerm$i.term$method) != 1) | !any(catTerm$i.term$method %in% i.meth.opt) ) } ){

# --> Make sure 'method' matches at least one of the possible method stopping rules.
    if( interactive() ){
      while( !( length(catTerm$i.term$method) == 1 & all(catTerm$i.term$method %in% i.meth.opt) ) )
        catTerm$i.term$method <- readline( paste("Select from ONLY ONE of the following info-based stopping rules - ",
    	                                         paste(i.meth.opt, collapse = ", "),
    	                                         ": ", sep = "")
    	                                 )    	                 	                    
    } else{
      stop( paste("'catTerm$i.term$method' must be ONLY ONE of ", paste(i.meth.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
    
  } # END if STATEMENT
 
 
## f2) i.term --> crit ##

## threshold ##
  if( any(catTerm$term == "info" & catTerm$i.term$method %in% c("threshold") ) ){
  
    if( length(catTerm$i.term$crit) != 1 )
      catTerm$i.term$crit <- NA
      
    if( any(is.na(catTerm$i.term$crit)) | any(!is.numeric(catTerm$i.term$crit)) ){
    	
# --> Make sure 'crit' is a positive number.
      if( interactive() ){
        repeat{
          catTerm$i.term$crit <- readline("Select a Fisher Information cut-off for variable termination: ")
          catTerm$i.term$crit <- suppressWarnings(as.numeric(catTerm$i.term$crit))
          
          if( !( all(!is.na(catTerm$i.term$crit)) & all(is.numeric(catTerm$i.term$crit)) ) )
            next;
          
          if( catTerm$i.term$crit >= 0 )
            break;
            
        } # END repeat LOOP
      } else{
        catTerm$i.term$crit <- 2
      } # END ifelse STATEMENTS
    } # END if STATEMENT

  } # END if STATEMENT

 
## change ##
  if( any(catTerm$term == "info" & catTerm$i.term$method %in% c("change") ) ){
  
    if( length(catTerm$i.term$crit) != 1 )
      catTerm$i.term$crit <- NA
      
    if( any(is.na(catTerm$i.term$crit)) | any(!is.numeric(catTerm$i.term$crit)) ){
    	
# --> Make sure 'crit' is a positive number.
      if( interactive() ){
        repeat{
          catTerm$i.term$crit <- readline("Select the maximum Fisher Information change for variable termination: ")
          catTerm$i.term$crit <- suppressWarnings(as.numeric(catTerm$i.term$crit))
          
          if( !( all(!is.na(catTerm$i.term$crit)) & all(is.numeric(catTerm$i.term$crit)) ) )
            next;
          
          if( catTerm$i.term$crit >= 0 )
            break;
            
        } # END repeat LOOP
      } else{
        catTerm$i.term$crit <- .5
      } # END ifelse STATEMENTS
    } # END if STATEMENT

  } # END if STATEMENT


   
## g) c.term ##

## g1) c.term --> method ## (only for the classification termination criterion)
  if( { any(catTerm$term == "class") &
  	    ( (length(catTerm$c.term$method) != 1) | !any(catTerm$c.term$method %in% c.meth.opt) ) } ){

# --> Make sure 'method' matches at least one of the possible classification stopping rules.   	    	    
    if( interactive() ){       			
      while( !( length(catTerm$c.term$method) == 1 & all(catTerm$c.term$method %in% c.meth.opt) ) )
        catTerm$c.term$method <- readline( paste("Select from ONLY ONE of the following classification-based stopping rules - ",
    	                                         paste(c.meth.opt, collapse = ", "),
    	                                         ": ", sep = "")
    	                                 )    	                 	                    
    } else{
      stop( paste("'catTerm$c.term$method' must be ONLY ONE of ", paste(c.meth.opt, collapse = ", "), sep = "" ) )
    } # END ifelse STATEMENTS
    
  } # END if STATEMENT


## g2) c.term --> bounds ## (only for the class termination criterion or weird selection)  
  if( any(catTerm$term == "class") | any(c(catStart$at, catMiddle$at) %in% "bounds") ){

# Step 1: Make sure that the classification bound(s) exist.
    if( !all(is.numeric(catTerm$c.term$bounds)) ){
  	
      if( interactive() ){      	
      	catTerm$c.term$bounds <- NA      	
      	while( !( all(is.numeric(catTerm$c.term$bounds)) & all(!is.na(catTerm$c.term$bounds)) ) ){
          cat("Select the classification bound(s) used for all simulees (press 'return' after each): ")
          catTerm$c.term$bounds <- suppressWarnings(as.numeric(scan(what = character(), quiet = TRUE)))
        } # END while LOOP        
      } else{
        catTerm$c.term$bounds <- 0
      } # END ifelse STATEMENT
      
    } # END if STATEMENT


# Step 2: Turn the classification bounds into a matrix.    
    if( !is.null( dim(catTerm$c.term$bounds) ) ){
    	
      if( (nrow(catTerm$c.term$bounds) != 1) & (nrow(catTerm$c.term$bounds) != dim(resp)[1]) )
        warning( "classification bounds are not specified for every simulee" )

# --> a) original matrix of bounds --> repeat that matrix for all people.      
      catTerm$c.term$bounds <- matrix(catTerm$c.term$bounds,
                                      nrow = dim(resp)[1],
                                      ncol = dim(catTerm$c.term$bounds)[2],
                                      byrow = TRUE)
        
    } else if( length(catTerm$c.term$bounds) == dim(resp)[1] ){
    	
# --> b) original vector of bounds equal to simulees --> one-column matrix.
      catTerm$c.term$bounds <- matrix(catTerm$c.term$bounds,
                                      nrow = dim(resp)[1],
                                      ncol = 1,
                                      byrow = TRUE)
    } else{
    	
# --> c) original vector of bounds NOT equal to simulees --> vector is on each row.
      catTerm$c.term$bounds <- matrix(catTerm$c.term$bounds,
                                      nrow = dim(resp)[1],
                                      ncol = length(catTerm$c.term$bounds),
                                      byrow = TRUE)
    } # END ifelse STATEMENTS
   
# --> d) sort the bounds (with a trick incase there is only one bound per person). 
    catTerm$c.term$bounds <- t( apply( cbind(catTerm$c.term$bounds, -Inf),
                                       MARGIN = 1, FUN = sort ) )[ , -1, drop = FALSE]

  } # END if STATEMENT
  
## g3) c.term --> categ ## (only for the class termination criterion)
  if( any(catTerm$term == "class") & is.null(catTerm$c.term$categ) )
    catTerm$c.term$categ <- 1:(ncol(catTerm$c.term$bounds) + 1)

  
## g4) c.term --> delta ## (only for the class termination criterion)
  if( { any( catTerm$term == "class" & catTerm$c.term$method %in% c("SPRT", "GLR") ) } ){
  	
    if( length(catTerm$c.term$delta) != 1 )
      catTerm$c.term$delta <- NA
    
    if( any(is.na(catTerm$c.term$delta)) | any(!is.numeric(catTerm$c.term$delta)) ){
 
# --> Make sure 'delta' is a positive number.    	
      if(interactive()){
        while( !( { length(catTerm$c.term$delta) == 1 &
        	        all(!is.na(catTerm$c.term$delta)) &
        	        all(is.numeric(catTerm$c.term$delta)) } ) ){
          catTerm$c.term$delta <- readline("Select an indifference region/testing half-width: ")
          catTerm$c.term$delta <- suppressWarnings(abs(as.numeric(catTerm$c.term$delta)))
        } # END while LOOP
      } else{
        catTerm$c.term$delta <- .1
      } # END ifelse STATEMENTS
       
    } # END if STATEMENT  
    
  } # END if STATEMENT


## g5) c.term --> alpha/beta ##
  if( any( catTerm$term == "class" & catTerm$c.term$method %in% c("SPRT", "GLR") ) ){
  
    if( length(catTerm$c.term$alpha) != 1 )
      catTerm$c.term$alpha <- NA
    if( length(catTerm$c.term$beta) != 1 )
      catTerm$c.term$beta <- NA
      
    if( any(is.na(catTerm$c.term$alpha)) | any(!is.numeric(catTerm$c.term$alpha)) ){

# --> Make sure 'alpha' is a positive number between 0 and 1.
      if( interactive() ){
        repeat{
          catTerm$c.term$alpha <- readline("Select an SPRT specified Type I error rate: ")
          catTerm$c.term$alpha <- suppressWarnings(as.numeric(catTerm$c.term$alpha))
          
          if( !( all(!is.na(catTerm$c.term$alpha)) & all(is.numeric(catTerm$c.term$alpha)) ) )
            next;
          
          if( catTerm$c.term$alpha < 1 & catTerm$c.term$alpha > 0 )
            break;
            
        } # END repeat LOOP
      } else{
        catTerm$c.term$alpha <- .05
      } # END ifelse STATEMENTS
    } # END if STATEMENTS
    
    if( any(is.na(catTerm$c.term$beta)) | any(!is.numeric(catTerm$c.term$beta)) ){
    	
# --> Make sure 'beta' is a positive number between 0 and 1.
      if( interactive() ){
        repeat{
          catTerm$c.term$beta <- readline("Select an SPRT specified Type II error rate: ")
          catTerm$c.term$beta <- suppressWarnings(as.numeric(catTerm$c.term$beta))
          
          if( !( all(!is.na(catTerm$c.term$beta)) & all(is.numeric(catTerm$c.term$beta)) ) )
            next;
          
          if( catTerm$c.term$beta < 1 & catTerm$c.term$beta > 0 )
            break;
            
        } # END repeat LOOP
      } else{
        catTerm$c.term$beta <- .05
      } # END ifelse STATEMENTS
    } # END if STATEMENT
    
  } # END if STATEMENT


## g6) c.term --> conf.lev ##
  if( any(catTerm$term == "class" & catTerm$c.term$method %in% c("CI") ) ){
  
    if( length(catTerm$c.term$conf.lev) != 1 )
      catTerm$c.term$conf.lev <- NA
      
    if( any(is.na(catTerm$c.term$conf.lev)) | any(!is.numeric(catTerm$c.term$conf.lev)) ){
    	
# --> Make sure 'conf.lev' is a positive number between 0 and 1.
      if( interactive() ){
        repeat{
          catTerm$c.term$conf.lev <- readline("Select a confidence interval level between 0 and 1: ")
          catTerm$c.term$conf.lev <- suppressWarnings(as.numeric(catTerm$c.term$conf.lev))
          
          if( !( all(!is.na(catTerm$c.term$conf.lev)) & all(is.numeric(catTerm$c.term$conf.lev)) ) )
            next;
          
          if( catTerm$c.term$conf.lev < 1 & catTerm$c.term$conf.lev > 0 )
            break;
            
        } # END repeat LOOP
      } else{
        catTerm$c.term$conf.lev <- .95
      } # END ifelse STATEMENTS
    } # END if STATEMENT

  } # END if STATEMENT


######################## END ARGUMENT CHECK SECTION #######################
###########################################################################


#~~~~~~~~~~~~~~~~~~~~~~~~#
# SETTING UP THE OBJECTS #
#~~~~~~~~~~~~~~~~~~~~~~~~#

# The number of people should be equal to the number of response matrix rows:   
  N <- dim(resp)[1]   # N is the number of row of theta/resp
  J <- dim(params)[1] # J is the number of rows of params

# A list to store individual, person attributes (final responses, final theta, etc.)
  cat_indiv  <- vector("list", length = N)

# Vectors to store CAT thetas, categories, test info, and SEM:
  cat_theta  <- vector("numeric",   length = N)
  cat_categ  <- vector("character", length = N)
  cat_info   <- vector("numeric",   length = N)
  cat_sem    <- vector("numeric",   length = N)
  cat_length <- vector("numeric",   length = N)
  cat_term   <- vector("character", length = N)
  
# Vectors to store TOT thetas, categories, test info, and SEM:
  tot_theta  <- vector("numeric",   length = N)
  tot_categ  <- vector("character", length = N)
  tot_info   <- vector("numeric",   length = N)
  tot_sem    <- vector("numeric",   length = N)
  
# Vectors to store TRUE thetas, categories:
  if( missing(theta) ){
    theta      <- NULL
    true_theta <- rep(NA, length = N)
  } else if( is.null(theta) ){
    true_theta <- rep(NA, length = N)
  } else{
    true_theta <- theta
  } # END ifelse STATEMENT
  
  true_categ <- vector("character", length = N)
  
  
# To store the selection rates (for Sympson-Hetter):
  S <- NULL


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# RUNNING THROUGH THE CAT: PER PERSON #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~#
# Progress Statement 1 #
#~~~~~~~~~~~~~~~~~~~~~~#

# For calibrations, we will want to suppress the upper progress bar.
  if( progress ){
  	cat("CAT Progress:\n")
    pb <- txtProgressBar(min = 0, max = N, initial = 0, char = "*", style = 3)
  }  # END if STATEMENT
    

# We will repeat the CAT for each person:
  for( i in 1:N ){
    
#####
# 1 # (INITIALIZING THE CAT)
#####

# First, the particular simulee's response vector:
    resp.i      <- resp[i, ]
    
    catStart.i  <- catStart
    catMiddle.i <- catMiddle
    catTerm.i   <- catTerm
  
# Second, the particular simulee's initial theta value and classification bounds:
    catStart.i$init.theta     <- catStart$init.theta[i]
    
    if( !is.null( dim(catTerm$c.term$bounds) ) )
      catTerm.i$c.term$bounds <- catTerm$c.term$bounds[i, ]

#####
# 2 # (THE FIRST SEVERAL ITEMS)
#####

# Responses, proximate theta ests/info/sem, item numbers/parameters of CAT:
    cat_resp.i  <- rep(NA, times = catTerm.i$n.max)
    cat_par.i   <- matrix(NA, nrow = catTerm.i$n.max, ncol = ncol(params))
    cat_it.i    <- cat_resp.i
  
    cat_theta.i <- c(catStart.i$init.theta, cat_resp.i)
    cat_info.i  <- cat_resp.i
    cat_sem.i   <- cat_resp.i
    
# Setting the classes, so the other functions will work.
    class(cat_par.i)  <- class(resp)
    class(cat_resp.i) <- mod
    
# A vector indicating whether or not the examinee has taken an item:
    it_flags    <- rep(0, nrow(params))
  
# Running the loop at the beginning of the CAT for each person:
    x <- startCat( params = params, resp = resp.i, mod = mod,
                   it_flags = it_flags,
                   catStart = catStart.i, catMiddle = catMiddle.i, catTerm = catTerm.i,
                   ddist = ddist, ... )
                 
# Updated to where we are in the CAT:
    j <- x$j
    S <- c(S, x$S)
  
#####
# 3 # (PRELIMINARY TERMINATION CHECK)
#####

# If the estimation method for early CATs was silly, then re-estimate.
    if( any(catStart.i$score %in% c("fixed", "step", "random")) ){
    	
# --> a) Build a character string for the WLE/EAP/BME estimation function,
      scoreFun <- paste(tolower(catMiddle.i$score), "Est", sep = "")
      
# --> b) Call the estimation function on stuff to this point,
      x  <- get(scoreFun)( resp = cat_resp.i[1:j],
                           params = cat_par.i[1:j, -c(1, ncol(params)), drop = FALSE],
                           range = catMiddle.i$range, mod = mod,
                           ddist = ddist, ... )
                           
      cat_theta.i[j + 1] <- x$theta
      cat_info.i[j]      <- x$info
      cat_sem.i[j]       <- x$sem
      
    } # END if STATEMENT

# Before entering the loop check the termination criterion:
    y <- termCat( params    = params[ , -c(1, ncol(params))],
                  resp      = resp.i,
                  mod       = mod,
                  it_flags  = it_flags,
                  cat_par   = cat_par.i[1:j , -c(1, ncol(cat_par.i)), drop = FALSE],
                  cat_resp  = cat_resp.i[1:j],
                  cat_theta = cat_theta.i[j + 1],
                  cat_info  = cat_info.i,
                  cat_sem   = cat_sem.i,
                  catStart  = catStart.i, catMiddle = catMiddle.i, catTerm = catTerm.i )
     
# And if we should stop, skip the while loop altogether:           
    stp <- y$cat_dec$stp
      
      
#####
# 4 # (THE REMAINING ITEMS)
#####

# We will repeat the next several steps until the CAT is over :(
    while( !stp ){
  	
      x <- middleCat( params    = params,
                      resp      = resp.i,
                      mod       = mod,
                      it_flags  = it_flags,
                      cat_par   = cat_par.i[1:j, , drop = FALSE],
                      cat_it    = cat_it.i[1:j],
                      cat_resp  = cat_resp.i[1:j],
                      cat_theta = cat_theta.i[j + 1],
                      cat_info  = cat_info.i,
                      cat_sem   = cat_sem.i,
                      catStart  = catStart.i, catMiddle = catMiddle.i, catTerm = catTerm.i,
                      ddist = ddist, ... )

# Updated to where we are in the CAT:
      j <- x$j
      S <- c(S, x$S)
      
      y <- termCat( params   = params[ , -c(1, ncol(params))],
                    resp     = resp.i,
                    mod      = mod,
                    it_flags = it_flags,
                    cat_par  = cat_par.i[1:j , -c(1, ncol(cat_par.i)), drop = FALSE],
                    cat_resp = cat_resp.i[1:j],
                    cat_theta = cat_theta.i[j + 1],
                    cat_info = cat_info.i,
                    cat_sem  = cat_sem.i,
                    catStart = catStart.i, catMiddle = catMiddle.i, catTerm = catTerm.i )
                    
      stp <- y$cat_dec$stp
      
    } # END while LOOP

#####
# 5 # (FINAL ASSIGNMENT AND FINAL CATEGORY STUFF)
#####

# Need to indicate: true  (a) class, (b) theta (if applicable), and
#                   total (a) class, (b) theta.

# First, find the appropriate theta estimate given resp and params.

# --> a) Build a character string for the WLE/EAP/BME estimation function,
    scoreFun <- paste(tolower(catTerm.i$score), "Est", sep = "")
      
# --> b) Call the estimation function on stuff to this point,
    x  <- get(scoreFun)( resp = resp.i,
                         params = params[ , -c(1, ncol(params))],
                         range = catMiddle.i$range, mod = mod, ddist = ddist, ... )
    
    tot_theta[i] <- x$theta
    tot_info[i]  <- x$info
    tot_sem[i]   <- x$sem
    
# If the EAP estimate didn't work, change estimation to BME with same dist/params:
    if( catTerm.i$score == "EAP" & is.nan(x$theta) ){

# --> a) Change scoring of the total theta to BME,
      catTerm$score <- "BME"
      
# --> b) Estimate the total thetas to this point as BME.
      x <- bmeEst( resp = resp[1:i, ],
                   params = params[ , -c(1, ncol(params))],
                   range = catMiddle.i$range, mod = mod, ddist = ddist, ... )
                   
      tot_theta[1:i] <- x$theta
      tot_info[1:i]  <- x$info
      tot_sem[1:i]   <- x$sem
      
      msg <- TRUE
    
    } # END if STATEMENT


# Second, determine the true and total classification.
  if( any(catTerm.i$term == "class") ){
  	
    tot_categ[i] <- catTerm.i$c.term$categ[sum( tot_theta[i] > catTerm.i$c.term$bounds ) + 1]  # based on total-test thetas
  	
    if( !is.null(theta) )
      true_categ[i] <- catTerm.i$c.term$categ[sum( theta[i] > catTerm.i$c.term$bounds ) + 1]  # based on true thetas
    
  } # END if STATEMENT


# Third, insert all of the individual CAT (thetas, sem, ... ) stuff
  cat_theta[i]        <- cat_theta.i[j + 1]
  cat_categ[i]        <- y$cat_categ
  cat_info[i]         <- cat_info.i[j]
  cat_sem[i]          <- cat_sem.i[j]
  cat_term[i]         <- y$cat_dec$term
  cat_length[i]       <- j
  
# Fourth, name the CAT parameters to put into the cat_indiv list.
  colnames(cat_par.i) <- colnames(params)
  
# Note that j had already been declared before another item/estimation.
  cat_indiv[[i]]      <- list( cat_theta  = cat_theta.i[1:(j + 1)],
                               cat_it     = cat_it.i[1:j],
                               cat_info   = c(NA, cat_info.i[1:j]),
                               cat_sem    = c(NA, cat_sem.i[1:j]),
                               cat_resp   = cat_resp.i[1:j],
                               cat_params = cat_par.i[1:j, ] )


#~~~~~~~~~~~~~~~~~~~~~~#
# Progress Statement 2 #
#~~~~~~~~~~~~~~~~~~~~~~#

# And indicate the progress thus far into the cat:
    if( progress )
      setTxtProgressBar(pb, value = i)

  } # END for i LOOP
  
  if( progress )
    cat("\n")


# And if we had to change EAP to BME, let people know:
  if( exists("msg", where = environment()) )
    if( msg )
      cat("\nToo large an item bank for full-test EAP estimation. Changed estimation for full-test to BME.\n")  
    
#~~~~~~~~~~~~~~~~~~~#
# FINISHING THE CAT #
#~~~~~~~~~~~~~~~~~~~#

# Determining the number of times each item was selected/administered:
  select <- tabulate(S, nbins = nrow(params))
  admin  <- tabulate(do.call(c, lapply(cat_indiv, FUN = function(x) x$cat_it)), nbins = nrow(params))
  names(select) <- names(admin) <- 1:nrow(params)

  ret <- list( cat_theta = cat_theta, cat_categ = cat_categ,      # CAT STUFF
               cat_info = cat_info, cat_sem = cat_sem,            # CAT STUFF
               cat_length = cat_length, cat_term = cat_term,      # CAT STUFF
               tot_theta = tot_theta, tot_categ = tot_categ,      # TOTAL TEST STUFF
               tot_info = tot_info, tot_sem = tot_sem,            # TOTAL TEST STUFF
               true_theta = true_theta, true_categ = true_categ,  # TRUE STUFF
               full_params = params, full_resp = resp,            # FULL PARAMETERS AND RESPONSES
               it_select = list(select = select, admin  = admin), # TABLE OF SELECTED AND ADMINISTERED ITEMS
               cat_indiv = cat_indiv,                             # INDIVIDUAL STUFF
               mod  = list(mod = mod, catStart = catStart, catMiddle = catMiddle, catTerm = catTerm) )

  class(ret) <- "catIrt"
  
  return(ret)
  
} # END catIrt FUNCTION