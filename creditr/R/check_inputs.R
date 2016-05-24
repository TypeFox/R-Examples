#' Check whether inputs from the data frame are valid.
#' 
#' \code{check_inputs} checks whether a data frame's inputs are valid. It is a 
#' minimum set of checks. Things such as recovery var are not checked, 
#' because some functions don't need them as input.
#' 
#' @inheritParams CS10
#'   
#' @return a data frame if not stopped by errors.
#'   
#' @examples 
#' x <- data.frame(date = as.Date(c("2014-04-22", "2014-04-22")),
#'                 currency = c("USD", "EUR"),
#'                 tenor = c(5, 5),
#'                 spread = c(120, 110),
#'                 coupon = c(100, 100),
#'                 recovery = c(0.4, 0.4),
#'                 notional = c(1000000, 1000000))
#' x <- check_inputs(x)

check_inputs <- function(x,
                         date.var     = "date", 
                         currency.var = "currency",
                         maturity.var = "maturity",
                         tenor.var    = "tenor",
                         spread.var   = "spread",
                         coupon.var   = "coupon",
                         notional.var = "notional",
                         notional     = 1000000,
                         recovery.var = "recovery",
                         recovery     = 0.4){

  stopifnot(is.data.frame(x))
  
  ## check if certain variables are contained in x
  
  stopifnot(c(currency.var, spread.var, coupon.var) %in% names(x))
  
  if(!notional.var %in% names(x)){
    x <- cbind(x, notional)
  }
  
  if(!recovery.var %in% names(x)){
    x <- cbind(x, recovery)
  }
  
  if(!date.var %in% names(x)){
    x <- cbind(x, date = Sys.Date())
  }
  
  if(!(tenor.var %in% names(x) | maturity.var %in% names(x))){
    x <- cbind(x, tenor = 5)
  }
  
  ## check if variables are defined in the correct classes
  
  ## Since R base does not have a is.date function for checking 
  ## whether a variable is of class "Date", we use
  ## inherits() instead. We don't want other packages.
  
  stopifnot(inherits(x[[date.var]], "Date"))
  
  x[[currency.var]] <- as.character(x[[currency.var]])
  
  ## check maturity OR tenor
  
  if(maturity.var %in% names(x) & tenor.var %in% names(x)){
    stop("do not provide both maturity and tenor")
  }else if(maturity.var %in% names(x)){
    stopifnot(inherits(x[[maturity.var]], "Date"))
  }else if(tenor.var %in% names(x)){
    stopifnot(is.numeric(x[[tenor.var]]))
  }else{
    stop("provide either tenor OR maturity")
  }
  
  stopifnot(is.numeric(x[[spread.var]])   & !(is.integer(x[[spread.var]])))
  stopifnot(is.numeric(x[[coupon.var]])   & !(is.integer(x[[coupon.var]])))
  stopifnot(is.numeric(x[[notional.var]]) & !(is.integer(x[[notional.var]])))
  
  ## change the column names anyway
  
  colnames(x)[which(colnames(x) == date.var)]     <- "date"
  colnames(x)[which(colnames(x) == currency.var)] <- "currency"
  colnames(x)[which(colnames(x) == maturity.var)] <- "maturity"
  colnames(x)[which(colnames(x) == tenor.var)]    <- "tenor"
  colnames(x)[which(colnames(x) == spread.var)]   <- "spread"
  colnames(x)[which(colnames(x) == coupon.var)]   <- "coupon"
  colnames(x)[which(colnames(x) == recovery.var)] <- "recovery"
  colnames(x)[which(colnames(x) == notional.var)] <- "notional" 

  return(x)
}
