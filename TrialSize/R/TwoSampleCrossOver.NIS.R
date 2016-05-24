TwoSampleCrossOver.NIS <-
function(alpha,beta,sigma,delta,margin){
n<-(qnorm(1-alpha)+qnorm(1-beta))^2*sigma^2/(2*(margin-delta)^2)
n
}
