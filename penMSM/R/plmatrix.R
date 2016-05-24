plmatrix <- function(psv, beta, constant){
  part1 <- dpenaltyfunction(psv=psv, beta=beta)
  part2 <- dapproxpenalty(psv=psv, beta=beta, constant=constant)
  ## ho <- as.numeric(psv%*%beta)
  ## if(ho == 0){
  ##     ho <- 1e-05}
  ## part2 <- part2/ho
  part3 <- psv%*%t(psv)
  result <- part1*part2*part3
  return(result)}
