WilliamsDesign.NIS <-
function(alpha,beta,sigma,sequence,delta,margin){
n<-(qnorm(1-alpha)+qnorm(1-beta))^2*sigma/(sequence*(delta-margin)^2)
n
}
