RelativeRiskCrossOver.NIS <-
function(alpha,beta,sigma,or,margin){
n<-(qnorm(1-alpha)+qnorm(1-beta))^2*sigma/(log(or)-margin)^2
n
}
