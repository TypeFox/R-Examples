#' Calculates Implied Recovery Rate
#' 
#' \code{implied_RR} that calculates the recovery rate implied by the CDS spread
#' and probability of default (pd) by using the ISDA model. This takes a data 
#' frame of inputs and returns a vector of the same length.
#' 
#' @inheritParams CS10
#' @param pd.var name of the column containing the probability of default rates.
#'   
#' @return implied recovery rate in percentage based on the general 
#'   approximation for a probability of default in the Bloomberg manual. The 
#'   actual calculation uses a complicated bootstrapping process, so the results
#'   may be marginally different.

implied_RR <- function(x, 
                       date.var     = "date",
                       tenor.var    = "tenor",
                       maturity.var = "maturity",
                       spread.var   = "spread",
                       pd.var       = "pd"){
  
  ## stop if the required columns are not contained in the data frame
  
  stopifnot(c(spread.var, pd.var, date.var) %in% names(x))
  
  ## stop if the required variables do not belong to the correct classes
  
  stopifnot(is.numeric(x[[spread.var]]))
  stopifnot(is.numeric(x[[pd.var]]))
  stopifnot(inherits(x[[date.var]], "Date"))
  
  if(maturity.var %in% names(x)){
  tenor <- as.numeric(x[[maturity.var]] - x[[date.var]])/360
  } else{
    if(tenor.var %in% names(x)){
      
      ## if we don't have maturity but we have tenor, we add a fake currency
      ## to x so we can use add_dates() to get the maturity date. Note the
      ## fake currency won't have any effect here, but just to pass the
      ## JPY check in add_dates
      
      x$currency <- "USD"
      x$maturity <- add_dates(x)$endDate
      x$currency <- NULL  ## remove the fake currency
      tenor <- x$tenor    ## just use the provided tenor by the user
    }else{
      stop("provide either tenor OR maturity")
    }
  }
  
  ## calculate the implied recovery rate 
  
  x <- c(100+((x[[spread.var]]*tenor/1e2)*(1/log(1-x[[pd.var]]))))
  
  return(x)
  
}
