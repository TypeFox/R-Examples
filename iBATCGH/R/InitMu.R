InitMu <-
function(deltak=c(-1,0,0.58,1),tauk=c(1,1,1,2),low_bounds=c(-Inf, -0.1, 0.1, 0.73) ,upp_bounds=c(-0.1, 0.1, 0.73, Inf)){
mu=c(rtnorm(1,mean=deltak[1],sd=tauk[1],lower=low_bounds[1],upper=upp_bounds[1]),rtnorm(1,mean=deltak[2],sd=tauk[2],lower=low_bounds[2],upper=upp_bounds[2]),rtnorm(1,mean=deltak[3],sd=tauk[3],lower=low_bounds[3],upper=upp_bounds[3]),rtnorm(1,mean=deltak[4],sd=tauk[4],lower=low_bounds[4],upper=upp_bounds[4]))
return (mu)
}
