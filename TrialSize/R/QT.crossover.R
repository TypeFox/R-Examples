QT.crossover <-
function(alpha, beta,pho,K,delta,gamma){
n=(pho+(1-pho)/K)*(qnorm(1-alpha/2)+qnorm(1-beta))^2/(delta^2-gamma*(qnorm(1-alpha/2)+qnorm(1-beta))^2)
}
