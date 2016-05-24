QT.PK.parallel <-
function(alpha, beta,pho,K,delta,v1,v2,tau1,tau2){
n=(2+(v1-v2)^2/(tau1^2+tau2^2))*(pho+(1-pho)/K)*(qnorm(1-alpha/2)+qnorm(1-beta))^2/delta^2
}
