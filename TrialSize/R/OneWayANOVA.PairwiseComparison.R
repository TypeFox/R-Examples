OneWayANOVA.PairwiseComparison <-
function(alpha,beta,tau,p1,p2,delta){
n.ij<-(qnorm(1-alpha/(2*tau))+qnorm(1-beta))^2*(p1*(1-p1)+p2*(1-p2))/delta^2
n.ij
}
