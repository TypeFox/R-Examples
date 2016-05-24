MeanWilliamsDesign.Equivalence <-
function(alpha,beta,sigma,k,delta,margin){
n<-(qnorm(1-alpha)+qnorm(1-beta/2))^2*sigma^2/(k*(delta-abs(margin))^2)
n
}
