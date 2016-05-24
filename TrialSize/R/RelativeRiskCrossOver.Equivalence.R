RelativeRiskCrossOver.Equivalence <-
function(alpha,beta,sigma,or,margin){
n<-(qnorm(1-alpha)+qnorm(1-beta/2))^2*sigma/(margin-abs(log(or)))^2
n
}
