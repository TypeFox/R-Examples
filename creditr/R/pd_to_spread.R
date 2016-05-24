#' Calculate spread with Default Probability
#' 
#' \code{pd_to_spread} to calculate spread using the probability of default,
#' tenor and recovery rate.
#' 
#' @inheritParams CS10
#' @param pd.var name of the column containing the probability of default in
#'   decimals.
#'   
#' @return vector containing the spread values in basis points, calculated by
#'   inverting the formula for probability of default given in the Bloomberg
#'   Manual
#'   
#' @seealso \code{\link{spread_to_pd}}

pd_to_spread <- function(x, 
                         recovery.var = "recovery", 
                         currency.var = "currency",
                         tenor.var    = "tenor", 
                         date.var     = "date", 
                         pd.var       = "pd"){
  
  stopifnot(is.data.frame(x))
    
  ## stop if the required variables are not contained in the data frame 
  
  stopifnot(c(recovery.var, tenor.var, currency.var, date.var, pd.var) 
            %in% names(x))
  
  ## stop if the variables do not belong to the correct class
  
  stopifnot(is.numeric(x[[tenor.var]]))
  stopifnot(is.numeric(x[[recovery.var]]))
  stopifnot(is.numeric(x[[pd.var]]))
  stopifnot(inherits(x[[date.var]], "Date"))
  
  x <- add_dates(x)
  
  ## calculate the exact time from the trade date till the maturity date Note:
  ## this 'time' is different from tenor. Let's say the trade date April 15,
  ## 2014 and the tenor is 5 years. Then the maturity date is June 20, 2019. So
  ## the time used to calculate the spread is the time between June 20, 2019 and
  ## April 15, 2014, which is 5.255556 years.

  time <- as.numeric(x$endDate - x$date) / 360
  
  ## calculate the spread by inverting the formula for probability of default
  
  spread <- 1e4*(x[[recovery.var]]-1)*log(1-x[[pd.var]])/time
  
  return(spread)
}
