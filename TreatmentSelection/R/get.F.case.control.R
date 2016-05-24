get.F.case.control <-
function(marker, event, trt, rho, return.fun = FALSE){

  #n <- length(marker)

  Y.D1 <- marker[event==1]
  Y.D0 <- marker[event==0]

  FY.D1 <- ecdf(Y.D1) #sum.I( Y, ">", Y.R1) #old way
  FY.D0 <- ecdf(Y.D0) #sum.I( Y, ">", Y.R0) #old way
  
  if(!return.fun){
  
    result <- FY.D1(marker)*rho[3] + FY.D0(marker)*(1-rho[3])
  
  }else {
 
    result <- function(x) FY.D1(x)*rho[3] + FY.D0(x)*(1-rho[3])

  }

  result

}
