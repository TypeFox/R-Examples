dapproxpenalty <- function(psv, beta, constant){
  currentpenalty <- as.numeric(t(psv)%*%beta)
  ## result <- currentpenalty/sqrt((currentpenalty^2)+constant)
  result <- 1/sqrt((currentpenalty^2)+constant)
  return(result)}
