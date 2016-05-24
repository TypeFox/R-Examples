ISCV.Equality <-
function(alpha,beta,CVt,CVr,m){

sigma1=1/(2*m)*CVt^2+CVt^4
sigma2=1/(2*m)*CVr^2+CVr^4

n=(sigma1+sigma2)(qnorm(1-alpha/2)+qnorm(1-beta))^2/(CVt-CVr)^2

}
