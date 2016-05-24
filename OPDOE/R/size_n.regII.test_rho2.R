size_n.regII.test_rho2 <- function(alpha=0.05,beta=0.2,delta=0.1){
  n <- ceiling(2*(qnorm(1-alpha/2)+qnorm(1-beta))^2/delta^2)+3
  n
}
