"simulate.gp" <-
function(object,...){
  # simulate sample processes from prior
  zero.coeff(object)  
  propose.coeff(object,block=0,proposal.sd=1)
  return(NULL)
}
