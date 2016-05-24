################################################################################
##
## $Id: backtest.function.R 1300 2008-08-27 21:01:11Z zhao $
##
## Returns an object of class backtest
##
################################################################################

## Returns a backtest object.  "x" is a data frame containing all raw
## data.  "in.var", "ret.var", "by.var", and "date.var" are character
## strings corresponding to column names in "x".  "in.var" and
## "ret.var" can be vectors containing multiple character strings.
## Only one "by.var" or "date.var" is allowed.  "date.var" currently
## acts the same as a "by.var".  "buckets" is a numeric specifying the
## number of buckets in which to group in.var and by.var,
## respectively.  "universe" is an expression for selecting a subset
## of "x".  "id.var" is only used in combination with "date.var"; it
## is a character string corresponding to a variable in "x" to be used
## to uniquely identify stocks across dates.

backtest <- function(x,
                     in.var,
                     ret.var,
                     universe,
                     by.var    = NULL,
                     date.var  = NULL,
                     id.var    = NULL,
                     buckets   = 5,
                     natural   = FALSE,
                     do.spread = TRUE,
                     by.period = TRUE,
                     overlaps  = 1){

  
  ## Corner Case: only one in.var of class factor allowed per backtest
  
  if(length(in.var) > 1){
    if(any(sapply(x[in.var], class) == "factor") ||
       any(sapply(x[in.var], class) == "character")){
        stop("Only one in.var of class factor or character allowed.")
    }
  }
 
  ## Corner Case: only one by.var allowed
  
  if(length(by.var) > 1){
    stop("Only one by.var allowed per backtest.")
  }

  ## Corner Case: only one id.var allowed

  if(length(id.var) > 1){
    stop("Only one id.var allowed per backtest.")
  }
  
  ## Corner Case: only one by.var allowed when using multiple ret.var

  if(length(ret.var) > 1 && (!is.null(by.var) || !is.null(date.var))){
    warning("Specifying by.var with multiple ret.vars is not supported. Proceed with caution.")
  }

  ## Must provide minimum of one in.var and one ret.var

  if(length(in.var) < 1 || length(ret.var) < 1){
    stop("At least one in.var and ret.var required.")
  }

  ## Only one date.var is allowed

  if(length(date.var) > 1){
    stop("Only one date.var is allowed.")
  }

  ## Natural backtests must have dates and ids.
  
  if(natural && (is.null(date.var) || is.null(id.var))){
    stop("Must specify date.var and id.var for a natural backtest.")
  }

  ## ret.var columns must be numeric
  
  if(!all(sapply(x[ret.var], class) == "numeric")){
    stop("All ret.var columns must be numeric")
  }
  
  ## Check "buckets"

  if(is.null(by.var)){
    buckets[2] <- 1
  }
  else{
    if(length(buckets) == 1){
      buckets[2] <- buckets[1]
    }
  }
  
  ## If overlaps is greater than 1, a date.var is required.

  if(overlaps > 1 && is.null(date.var)){
    stop("If overlaps is greater than 1, a date.var is required.")
  }

  ## If overlaps is greater than 0, in.var must be of length one and numeric.

  if(overlaps > 1 && (length(in.var) != 1 || !is.numeric(x[[in.var]]))){
    stop("If overlaps is greater than 1, in.var must be of length one and numeric.")
  }

  ## Overlap option must include id.var
  
  if(overlaps > 1 && length(id.var) < 1){
    stop("If overlaps is greater than 1, there must be an id.var.")
  }


  ## Overlaps must be less than the number of periods

  if(overlaps != 1 && overlaps > length(unique(x[[date.var]]))){
    stop("Overlaps must be less than the number of periods.")
  }

  ## If overlaps > 1, can only use one ret.var
  
  if(overlaps > 1 && length(ret.var) > 1){
    stop("The multiple overlap option can only accept one ret.var.")
  }
  
  ## At a minimum, the length of in.var must be greater than the
  ## number of buckets.  Otherwise, we will get an error when we try
  ## to create more quantiles then there are observations.

  if(any(sapply(x[in.var], function(x){sum(!is.na(x))}) < (buckets[1] - 5))){
    stop("The number of non-NA in.var must be at least 5 more than the number of buckets.") 
  }
  
  ## The length of by.var must also be greater than the number of
  ## buckets to properly create quantiles
  
  if(!is.null(by.var) && sum(!is.na(x[by.var])) < (buckets[2] - 5)){
    stop("The number of non-NA by.var must be at least 5 more than the number of buckets.")
  }

  ## Check for by.period and date

  if(isTRUE(by.period) && is.null(date.var)){
    stop("date.var required if by.period = TRUE (the default)")
  }

  
  ## Save by.var for by.var slot
  
  by.specified <- by.var

  ## Use by.var or date.var?
  
  if(!is.null(date.var)){

    if(!is.null(by.var)){
      stop("Cannot specify both by.var and date.var.")
    }
    else{
      by.var <- date.var
    }
  }
  
  ## Evaluate "universe"
  
  if(!missing(universe)){
    univ <- eval(substitute(universe), x, parent.frame())
    univ <- univ & !is.na(univ)
    x    <- x[univ,]
  }
  
  ## Attach "by.var" factor to "x"
  
  if(is.null(by.var)){
    x$by.factor <- rep(factor(1), times=nrow(x))
  }
  else{
    if(is.null(date.var)){
      x$by.factor <- categorize(x[[by.var]], n = buckets[2])
    }
    else{
      x$by.factor <- categorize(x[[by.var]], n = buckets[2], is.date = TRUE)
    }
  }
  
  ## Factor a character in.var
  
  if(length(in.var) == 1 && is.character(x[in.var]))
    x[in.var] <- as.factor(x[in.var])
  
  ## Make sure in.var is a factor / factorized... by looping over in.var
  
  in.factor <- data.frame(array(dim = c(nrow(x), length(in.var)), dimnames = list(NULL,in.var)))
  
  ## Intialize weights
  
  x$weight <- 1
  
  ## Build buckets
  
  for(i in in.var){
    
    ## If in.var is not a factor (is numeric)

    if(!is.factor(x[[i]])){
      
      if(by.period && !all(tapply(x[[i]], x[[date.var]],
                                          function(x){
                                            if(sum(!is.na(x)) < buckets[1]){
                                              return(FALSE)
                                            }
                                            return(TRUE)
                                          }
                                          )
                                   )
         ){
        stop("Not enough observations to fill each bucket by period.")
      }
      
      ## If we are doing by period, form buckets for every date 

      if(by.period){
        in.factor[[i]] <- as.factor(unsplit(lapply(split(x[[i]], x[[date.var]]),
                                                   function(x){
                                                     categorize(x, n = buckets[1])
                                                   }),
                                            x[[date.var]]))
      }
      
      ## If we are not doing by period bucketizing, form all buckets simultaneously

      else{
        in.factor[[i]] <- categorize(x[[i]], n = buckets[1])
      }
      
      
      ## Name levels, the lowest bucket is the short part of the
      ## portfolio, the highest bucket is the long part      

      levels(in.factor[[i]])[1]          <- "low"
      levels(in.factor[[i]])[buckets[1]] <- "high"
      
      ## Make sure there is data in every bucket
      
      if(length(levels(in.factor[[i]])) != buckets[1]){
        stop(paste("Encountered quantiles with no observations.  This can",
                   "occur with very little data or very regular",
                   "(usually synthesized) data."))
      }
          
      ## Recalculate weights based on overlaps

      if(overlaps > 1){
        levels(in.factor[[i]])[1]                  <- "low"
        levels(in.factor[[i]])[buckets[1]]         <- "high"
        levels(in.factor[[i]])[2:(buckets[1] - 1)] <- "mid"

        in.factor.col <- paste("in.factor.", i, sep = "")
        
        x[[in.factor.col]] <- in.factor[[i]]
        
        x <- overlaps.compute(x, in.factor.col, date.var, id.var, overlaps)
      }
    }

    
    ## If the in.var is a factor, just set in.factor to the in.var

    else{
      in.factor[[i]] <- x[[i]]
    }
  }
  

  
  invisible(backtest.compute(x, in.factor, ret.var, by.var, date.var,
                             natural, by.specified, do.spread, id.var,
                             by.period, overlaps)) 
}
