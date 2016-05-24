TwoSampleSeqCrossOver.Equality <-
function(alpha,beta,sigma,sequence,delta){
n<-(qnorm(1-alpha/2)+qnorm(1-beta))^2*sigma/(sequence*delta^2)
n
}
