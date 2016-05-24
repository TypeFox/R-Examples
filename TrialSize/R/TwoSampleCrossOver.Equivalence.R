TwoSampleCrossOver.Equivalence <-
function(alpha,beta,sigma,delta,margin){
if (delta==0) n<-(qnorm(1-alpha)+qnorm(1-beta/2))^2*sigma^2/(2*delta^2)
else n<-(qnorm(1-alpha)+qnorm(1-beta))^2*sigma^2/(2*(delta-abs(margin))^2)
n
}
