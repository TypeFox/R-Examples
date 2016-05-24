Dose.Min.Effect <-
function(alpha, beta,qt,sigma,delta){

#  not real t distribution qt(1-alpha,K) is from the table 12.1.1-12.1.4 in the book 
n=2*sigma^2*(qt+qnorm(1-beta))^2/delta^2
}
