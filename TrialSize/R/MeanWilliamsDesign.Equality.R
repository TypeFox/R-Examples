MeanWilliamsDesign.Equality <-
function(alpha,beta,sigma,k,margin){
n<-(qnorm(1-alpha/2)+qnorm(1-beta))^2*sigma^2/(k*margin^2)
n
}
