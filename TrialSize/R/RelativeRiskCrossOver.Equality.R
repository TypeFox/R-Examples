RelativeRiskCrossOver.Equality <-
function(alpha,beta,sigma,or){
n<-(qnorm(1-alpha/2)+qnorm(1-beta))^2*sigma/(log(or))^2
n
}
