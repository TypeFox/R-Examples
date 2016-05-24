PBE <-
function(alpha, beta, sigma1.1,sigmatt,sigmatr,sigmabt,sigmabr,rho,a,delta,lamda){

temp=2*delta^2*sigma1.1^2+sigmatt^4+(1+a)^2*sigmatr^4-2*(1+a)*rho^2*sigmabt^2*sigmabr^2
n=temp*(qnorm(1-alpha)+qnorm(1-beta))^2/(lamda^2)
return(n)
}
