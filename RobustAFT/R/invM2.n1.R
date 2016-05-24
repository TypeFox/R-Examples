"invM2.n1" <-
function(u,theta,sigma,rs,wi,estim=c("SA","TMLA","SI","TMLI")) {
p <- 1; n <- length(rs); xk <- 1.5477
if (estim=="SA")  {zpsp <- psp.weight(rs,ips=2,xk=xk)
                   a1c  <- sum(zpsp)/(n*sigma)
                   b1c  <- sum(zpsp*rs)/(n*sigma)
                   zpsi <- psi.weight(rs,ips=2,xk=xk)
                   a2c  <- sum(zpsi)/(n*sigma)
                   b2c  <- sum(zpsi*rs)/(n*sigma)}
if (estim=="TMLA"){a1c  <- sum(wi)/(n*sigma)
                   b1c  <- sum(rs*wi)/(n*sigma)
                   a2c  <- 2*b1c
                   b2c  <- 2*sum(rs^2*wi)/(n*sigma)}
if (estim=="SI")  {a1c   <- integrate(Pspphi.n, lower=-xk,upper=xk)$value/sigma
                   b2c   <- integrate(Psizphi.n,lower=-xk,upper=xk)$value/sigma
                   a2c   <- 0; b1c <- 0}
if (estim=="TMLI"){a1c   <- (2*pnorm(u)-1)/sigma
                   b2c   <- 4*(-u*dnorm(u)+pnorm(u)-0.5)/sigma
                   a2c   <- 0; b1c <- 0}
M  <- matrix(0,ncol=p+1,nrow=p+1)
M[1,1] <- a1c; M[1,2] <- b1c
M[2,1] <- a2c; M[2,2] <- b2c
Minv <- solve(M)
list(Minv=Minv)}

