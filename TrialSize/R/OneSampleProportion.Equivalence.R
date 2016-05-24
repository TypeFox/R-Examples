OneSampleProportion.Equivalence <-
function(alpha,beta,p,delta,margin){
n<-(qnorm(1-alpha)+qnorm(1-beta/2))^2*p*(1-p)/(margin-abs(delta))^2
n
}
