Dose.Response.Linear <-
function(alpha, beta, sigma,mui,ci,fi){
epsilon=sum(ci*mui)
n=(qnorm(1-alpha)+qnorm(1-beta))^2*sigma^2*sum(ci^2/fi)/(epsilon^2)
}
