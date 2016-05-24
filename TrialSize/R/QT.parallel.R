QT.parallel <-
function(alpha, beta,pho,K,delta){
n=2*(pho+(1-pho)/K)*(qnorm(1-alpha/2)+qnorm(1-beta))^2/delta^2
}
