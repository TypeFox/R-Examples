QT.PK.crossover <-
function(alpha, beta,pho,K,delta,gamma,v1,v2,tau1,tau2){
n=(1+(v1-v2)^2/(tau1^2+tau2^2))*(pho+(1-pho)/K)*(qnorm(1-alpha/2)+qnorm(1-beta))^2/(delta^2-gamma*(qnorm(1-alpha/2)+qnorm(1-beta))^2)
}
