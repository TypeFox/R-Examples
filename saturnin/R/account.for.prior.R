account.for.prior <-
function(prob,q0){
  p0 <- 2/nrow(prob)
  ratio <- (p0*(1-q0))/(q0*(1-p0))
  prob.q0 <- apply(prob,
                    1:2,
                    function(y) (1 + ratio*((1-y)/y))^{-1})
  return(prob.q0)
}
