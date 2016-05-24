OneSampleProportion.Equality <-
function(alpha,beta,p,delta){
n<-(qnorm(1-alpha/2)+qnorm(1-beta))^2*p*(1-p)/delta^2
n
}
