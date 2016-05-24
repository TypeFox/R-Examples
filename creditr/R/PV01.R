#' Calculate PV01
#' 
#' \code{PV01} to calculate present value 01 or present value of a stream of 1bp
#' payments
#' 
#' @name PV01
#'   
#' @inheritParams CS10
#' @param principal.var name of the column containing the principal or clean
#'   upfront values of the CDS
#'   
#' @return Vector containing the PV01 values

PV01 <- function(x, 
                 principal.var = "principal", 
                 spread.var = "spread", 
                 coupon.var = "coupon", 
                 notional.var = "notional"){
  
  stopifnot(c(principal.var, spread.var, coupon.var, notional.var) %in% names(x))
  
  stopifnot(is.numeric(x[[principal.var]]))
  stopifnot(is.numeric(x[[notional.var]]))
  stopifnot(is.numeric(x[[spread.var]]))
  stopifnot(is.numeric(x[[coupon.var]]))
  
  PV01 <- (abs(x[[principal.var]]) / x[[notional.var]])*
    (10000 / abs(x[[spread.var]] - x[[coupon.var]]))
  
  return(PV01)
}