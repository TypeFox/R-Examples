WilliamsDesign.Equivalence <-
function(alpha,beta,sigma,sequence,delta,margin){
n<-(qnorm(1-alpha)+qnorm(1-beta/2))^2*sigma/(sequence*(margin-abs(delta))^2)
n
}
