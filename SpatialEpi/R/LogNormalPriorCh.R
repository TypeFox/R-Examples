LogNormalPriorCh <-
function(theta1, theta2, prob1, prob2){	
zq1 <- qnorm(prob1)
zq2 <- qnorm(prob2)
mu <- log(theta1)*zq2/(zq2-zq1) - log(theta2)*zq1/(zq2-zq1)
sigma <- (log(theta1)-log(theta2))/(zq1-zq2)
list(mu=mu,sigma=sigma)
}
