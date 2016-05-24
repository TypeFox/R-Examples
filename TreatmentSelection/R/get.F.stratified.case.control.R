get.F.stratified.case.control <-
function(marker, event, trt, rho, return.fun = FALSE){
 
  #n <- length(marker)     

  #rho[1] = Pr( event = 1 | trt = 0 )
  #rho[2] = Pr( event = 1 | trt = 1 )

  Pr.D1.trt1 <- rho[5]
  Pr.D0.trt1 <- rho[4]
  Pr.D1.trt0 <- rho[3]
  Pr.D0.trt0 <- rho[2]
  
  Y.D1.trt1 <- marker[event==1 & trt==1]
  Y.D0.trt1 <- marker[event==0 & trt==1]
  Y.D1.trt0 <- marker[event==1 & trt==0]
  Y.D0.trt0 <- marker[event==0 & trt==0]

  
                                  #old way
  FY.D1.trt1 <- ecdf(Y.D1.trt1) #sum.I( Y, ">", Y.R1.trt1)
  FY.D0.trt1 <- ecdf(Y.D0.trt1) #sum.I( Y, ">", Y.R0.trt1)
  FY.D1.trt0 <- ecdf(Y.D1.trt0) #sum.I( Y, ">", Y.R1.trt0)
  FY.D0.trt0 <- ecdf(Y.D0.trt0) #sum.I( Y, ">", Y.R0.trt0)    


  if(!return.fun){
  
    result <- FY.D1.trt1(marker)*(Pr.D1.trt1) + FY.D0.trt1(marker)*(Pr.D0.trt1) + FY.D1.trt0(marker)*(Pr.D1.trt0) + FY.D0.trt0(marker)*(Pr.D0.trt0) 

  }else{

      result <- function(x) FY.D1.trt1(x)*(Pr.D1.trt1) + FY.D0.trt1(x)*(Pr.D0.trt1) + FY.D1.trt0(x)*(Pr.D1.trt0) + FY.D0.trt0(x)*(Pr.D0.trt0) 
  
  }


  result
}
