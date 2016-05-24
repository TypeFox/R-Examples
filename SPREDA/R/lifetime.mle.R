lifetime.mle <-
function(dat, minusloglik, starts, method = "BFGS",hessian = TRUE,...)
{
  call=match.call()
  f = function(p) {
    minusloglik(dat,p) 
    }
  oout = optim(starts, f, method = method, hessian = hessian,...)#,control=list(trace=T))
  coef = oout$par
  #browser()
  if(hessian)
  {
   vcov =solve(oout$hessian)
  }else{
         vcov=NULL
       }
  min = oout$value
  invisible(list(call = call, coef = coef,vcov = vcov, min = min,dat=dat,minusloglik=minusloglik))
}
