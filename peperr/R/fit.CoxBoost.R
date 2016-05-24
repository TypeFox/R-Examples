`fit.CoxBoost` <-
function(response, x, cplx, ...){
   require(CoxBoost)
   time <- response[,"time"]
   status <- response[,"status"]

#    if (is.list(cplx)){
#    step <- cplx$stepno
#    pen <- cplx$stepsize.factor
#    res <- CoxBoost(time, status, x, stepno=step, stepsize.factor=pen, ...)
#    } else {
   res <- CoxBoost(time, status, x, stepno=cplx, ...)
#   }
   res
}

