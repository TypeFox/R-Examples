work.kriging <-
function(Gamma,gamma,z)
{
 ## Gamma : ordinary kriging matrix
 ## gamma : semivariance between query point and neighbors
 ## z     : observed values
 G <- rbind(cbind(Gamma,1),c(rep(1,ncol(Gamma)),0))
 g <- c(gamma,1)
 lambda <- solve(G,g)
 tmp <- lambda[-nrow(G)]
 m <- lambda[nrow(G)]
 zhat <- sum(tmp * z)
 sigmasq <- m + sum(tmp * gamma)
 sigmasq <- pmax(0,sigmasq)  # ordinary kriging may be negative
 out <- c(zhat,sqrt(sigmasq))
 names(out) <- c('fit','se')
 out
}
