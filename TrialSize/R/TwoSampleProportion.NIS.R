TwoSampleProportion.NIS <-
function(alpha,beta,p1,p2,k,delta,margin){
n2<-(qnorm(1-alpha)+qnorm(1-beta))^2*(p1*(1-p1)/k+p2*(1-p2))/(delta-margin)^2
n1<-k*n2
n1
}
