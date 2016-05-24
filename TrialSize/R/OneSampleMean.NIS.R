OneSampleMean.NIS <-
function(alpha,beta,sigma,margin,delta){
n<-(qnorm(1-alpha)+qnorm(1-beta))^2*sigma^2/(margin-delta)^2
n
}
