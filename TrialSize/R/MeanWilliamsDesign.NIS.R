MeanWilliamsDesign.NIS <-
function(alpha,beta,sigma,k,delta,margin){
n<-(qnorm(1-alpha)+qnorm(1-beta))^2*sigma^2/(k*(margin-delta)^2)
n
}
