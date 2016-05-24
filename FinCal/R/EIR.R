#' Equivalent/proportional Interest Rates
#' @description An interest rate to be applied n times p.a. can be converted to an equivalent rate to be applied p times p.a.
#' @param r interest rate to be applied n times per year (r is annual rate!)
#' @param n times that the interest rate r were compounded per year
#' @param p times that the equivalent rate were compounded per year
#' @param type equivalent interest rates ('e',default) or proportional interest rates ('p')
#' @export
#' @examples
#' # monthly interest rat equivalent to 5% compounded per year
#' EIR(r=0.05,n=1,p=12)
#' 
#' # monthly interest rat equivalent to 5% compounded per half year
#' EIR(r=0.05,n=2,p=12)
#' 
#' # monthly interest rat equivalent to 5% compounded per quarter
#' EIR(r=0.05,n=4,p=12)
#' 
#' # annual interest rate equivalent to 5% compounded per month
#' EIR(r=0.05,n=12,p=1)
#' # this is equivalent to
#' ear(r=0.05,m=12)
#' 
#' # quarter interest rate equivalent to 5% compounded per year
#' EIR(r=0.05,n=1,p=4)
#' 
#' # quarter interest rate equivalent to 5% compounded per month
#' EIR(r=0.05,n=12,p=4)
#' 
#' # monthly proportional interest rate which is equivalent to a simple annual interest
#' EIR(r=0.05,p=12,type='p')
EIR <- function(r,n=1,p=12,type=c("e", "p")){
  type = match.arg(type)
  if(type == "e"){
    eir=(1+r/n)^(n/p)-1
  }else if(type == "p"){
    eir=r/p
  }else{
    stop("type must be 'e' or 'p'")
  }
  return(eir)
}