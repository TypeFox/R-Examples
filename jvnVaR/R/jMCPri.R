jMCPri <-
function(s0, mu,sigma,m){
# h is the step of sim
if (m==0){
return(s0)
}
n <- 100
d <- 1/n
mat <- matrix(rep(1,m*n), nrow = m, ncol = n)
dri <- d * mat
vol <- sqrt(abs(d)) * mat
z <- matrix(rnorm(m*n), nrow = m, ncol = n)
s <- 1 + mu * dri + sigma * vol * z
s[s<0] <- 1
s <- s0 * apply(s,1,prod)
return(s)
}
