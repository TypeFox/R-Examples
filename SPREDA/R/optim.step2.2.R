optim.step2.2 <-
function(Dt,  A, B, r1, r2){
  num=-A*exp(r1)
  den=1+exp(-Dt/(B*exp(r2)))
  res=num/den
  return(res)
}
