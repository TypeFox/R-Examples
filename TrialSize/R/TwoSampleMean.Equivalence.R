TwoSampleMean.Equivalence <-
function(alpha,beta,sigma,k,delta,margin){
n2<-(qnorm(1-alpha)+qnorm(1-beta/2))^2*sigma^2*(1+1/k)/(delta-abs(margin))^2
n2
n1<-k*n2
}
