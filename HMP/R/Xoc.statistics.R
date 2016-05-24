Xoc.statistics <-
function(group.data, epsilon=10^(-4)){
fit <- lapply(group.data, function(x, epsilon){
dirmult::dirmult(x, initscalar=DM.MoM(x)$theta, epsilon=epsilon, trace=FALSE)
}, epsilon=epsilon)

thetas <- unlist(lapply(fit, function(x){x$theta}))
inipi <- lapply(fit, function(x) x$pi)

logliks <- unlist(lapply(fit, function(x){x$loglik}))
fit1 <- dirmult::equalTheta(group.data, theta=mean(thetas), epsilon=epsilon, trace=FALSE, initPi=inipi)
Xoc <- as.vector(-2*(fit1$loglik-sum(logliks)))

return(Xoc)
}
