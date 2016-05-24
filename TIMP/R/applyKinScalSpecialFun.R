"applyKinScalSpecialFun" <- function(ind, theta, parlist, spec) {
  # spec should contain an argument "functionlist"
  # and an additional list "inputs"

  inputs <- spec$inputs 
  funlist <- spec$functionlist 

  parvec <- parlist[[ind]]
  funx <- funlist[[ind]]
  inputsfixed <- inputs[[ind]] 

  kall <- theta 
  res <- funx(kall, parvec, inputsfixed)
  res
  
}
