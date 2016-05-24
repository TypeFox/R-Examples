Nonpara.Two.Sample <-
function(alpha, beta,k, p1,p2,p3){
n=(qnorm(1-alpha/2)*sqrt(k*(k+1)/12)+qnorm(1-beta)*sqrt(k^2*(p2-p1^2)+k*(p3-p1^2)))^2/(k^2*(1/2-p1)^2)
}
