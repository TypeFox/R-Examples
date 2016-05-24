OneSampleProportion.NIS <-
function(alpha,beta,p,delta,margin){
n<-(qnorm(1-alpha)+qnorm(1-beta))^2*p*(1-p)/(delta-margin)^2
n
}
