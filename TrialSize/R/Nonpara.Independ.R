Nonpara.Independ <-
function(alpha, beta,p1,p2){
n=4*(qnorm(1-alpha/2)/3+qnorm(1-beta)*sqrt(2*p2-1-(2*p1-1)^2))^2/(2*p1-1)^2
}
