#' Calcualte Default Probability with Spread
#' 
#' \code{spread_to_pd} approximates the default probability at time given the 
#' spread
#' 
#' @inheritParams CS10
#'   
#' @return vector containing the probability of default, calculated by using the
#'   formula for probability of default given in the Bloomberg Manual
#'   
#' @seealso \code{\link{pd_to_spread}}

spread_to_pd <- function(x, 
                         recovery.var = "recovery", 
                         currency.var = "currency",
                         tenor.var    = "tenor",
                         maturity.var = "maturity",
                         date.var     = "date", 
                         spread.var   = "spread"){
  
  stopifnot(is.data.frame(x))
  
  has.maturity.var <- !(is.null(maturity.var) || is.null(x[[maturity.var]]))
  has.tenor.var    <- !(is.null(tenor.var)    || is.null(x[[tenor.var]]))
  stopifnot(xor(has.maturity.var, has.tenor.var))
  
  ## stop if the required variables are not contained in the data frame 
  
  stopifnot(c(recovery.var, tenor.var, currency.var, date.var, spread.var) 
            %in% names(x))
  
  ## Coerce to desired classes
  
  x[[currency.var]] <- as.character(x[[currency.var]])
  x[[date.var]] <- as.Date(x[[date.var]])
  
  ## stop if the variables do not belong to the correct class
  
  stopifnot(!has.tenor.var || is.numeric(x[[tenor.var]]))
  stopifnot(is.numeric(x[[recovery.var]]))
  stopifnot(is.numeric(x[[spread.var]]))
  stopifnot(inherits(x[[date.var]], "Date"))
  
  x <- add_dates(x,
                 date.var     = date.var,
                 tenor.var    = tenor.var,
                 maturity.var = maturity.var,
                 currency.var = currency.var)
  
  ## calculate the exact time from the trade date till the maturity date Note:
  ## this 'time' is different from tenor. Let's say the trade date April 15,
  ## 2014 and the tenor is 5 years. Then the maturity date is June 20, 2019. So
  ## the time used to calculate the spread is the time between June 20, 2019 and
  ## April 15, 2014, which is 5.255556 years.
  
  time <- as.numeric(x$endDate - x$date)/360
  
  ## Bloomberg Approximation
  
  pd <- 1 - exp(-x[[spread.var]] / 1e4 * time / (1 - x[[recovery.var]]))
  
  return(pd)
}