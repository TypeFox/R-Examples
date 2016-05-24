jMCPriLim <-
function(s0,L,U, mu,sigma,m){
# h is the step of sim
if (m==0){
return(s0)
}
n <- 100
d <- 1/n
dri <- d * rep(1,m)
vol <- sqrt(abs(d)) * rep(1,m)
z <- matrix(rnorm(m*n), nrow = m, ncol = n)
s <- rep(1,m)
for (i in 1:n){ 
s <- s*(1 + mu * dri + sigma * vol * z[,i])
s[s<1+L] <- 1+L
s[s>1+U] <- 1+U
}
return(s0*s)
}
