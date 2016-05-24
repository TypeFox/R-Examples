Nonpara.One.Sample <-
function(alpha, beta, p2,p3,p4){
n=(qnorm(1-alpha/2)/sqrt(12)+qnorm(1-beta)*sqrt(p3+4*p4-4*p2^2))^2/(1/4-p2)^2
}
