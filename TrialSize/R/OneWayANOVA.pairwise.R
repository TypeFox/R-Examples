OneWayANOVA.pairwise <-
function(alpha,beta,tau,sigma,margin
){
n.ij<-2*(qnorm(1-alpha/(2*tau))+qnorm(1-beta))^2*sigma^2/margin^2
n.ij
}
