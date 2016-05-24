ISCV.NIS <-
function(alpha,beta,CVt,CVr,m,margin){

sigma1=1/(2*m)*CVt^2+CVt^4
sigma2=1/(2*m)*CVr^2+CVr^4

n=(sigma1+sigma2)*(qnorm(1-alpha)+qnorm(1-beta))^2/(CVt-CVr-margin)^2

}
