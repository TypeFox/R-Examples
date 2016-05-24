"ldmvnorm" <-
function(Y,S) {     # log density of a matrix with mvn rows
 n<-dim(Y)[1]
 p<-dim(Y)[2]
 -.5*n*log(det(S)) -.5*n*p*log(2*pi)-.5*sum( diag( solve(S)%*%t(Y)%*%Y)) }

