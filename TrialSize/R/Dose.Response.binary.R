Dose.Response.binary <-
function(alpha, beta,pi,ci,fi){
epsilon=sum(ci*pi)
pmean=mean(pi)
n=(qnorm(1-alpha)*sqrt(sum(ci^2*pmean*(1-pmean)/fi))+qnorm(1-beta)*sqrt(sum(ci^2*pi*(1-pi)/fi)))^2/(epsilon^2)
}
