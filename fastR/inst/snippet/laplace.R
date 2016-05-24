dlaplace <- function(x,theta,lambda) {
    0.5 * lambda * exp(-lambda * abs(x-theta))
}
# two ways to do plaplace:
integrate(function(x) {dlaplace(x,1,2)},-Inf,Inf)     # should = 1
plaplace1 <- function(q,theta=0,lambda=1) {
    integrate(function(x) {dlaplace(x,theta,lambda)},-Inf,q)$value
}

plaplace2 <- function(q,theta=0,lambda=1) {
    if (q < theta) return( 0.5 * (1-pexp(theta-q,lambda)) )
    return( 0.5 + 0.5 * pexp(q-theta,lambda) )
    }
# should get same results either way:
plaplace1(3,2,1)
plaplace1(3,2,1) - plaplace1(-3,2,1)
plaplace2(3,2,1)
plaplace2(3,2,1) - plaplace2(-3,2,1)
x <- c(1.00,-1.43,0.62,0.87,-0.66,-0.59,1.30,-1.23,-1.53,-1.94)
loglik1 <- function(theta, x) {
    m <- theta[1]; lambda <- theta[2]
    return( sum( log(0.5) + dexp(abs(x-m),rate=lambda, log=T)) )
}
mle <- nlmax(loglik1,p=c(0,1),x=x)$estimate; mle

# method of moments estimates
# estimate for theta is the sample mean since E(X) == est.theta:
est.theta = mean(x); est.theta
# estimate for variance satisfies v == 2/est.lambda^2:
v = var(x) * (n-1) / n
est.lambda = sqrt(2 / v); est.lambda
