gauss.copula <-
function(rho,U,V){
  qU <- qnorm(U,0,1)
  qV <- qnorm(V,0,1)
  exp((-rho^2*(sum(qU^2)+sum(qV^2))+2*rho*sum(qU*qV))/(2*(1-rho^2)))/(sqrt(1-rho^2))
}
