TwoSampleProportion.Equality <-
function(alpha,beta,p1,p2,k,delta){
n2<-(qnorm(1-alpha/2)+qnorm(1-beta))^2*(p1*(1-p1)/k+p2*(1-p2))/delta^2
n1<-k*n2
n1
}
