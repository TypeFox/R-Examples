"invM2.w1" <-
function(l,u,theta,sigma,rs,wi,estim=c("SA","TMLA","SI","TMLI")) {
n <- length(rs); xk <- 1.717817
if (estim=="SA")  {zpsp <- psp.weight(rs,ips=2,xk=1.717817)
                   a1c  <- sum(zpsp)/(n*sigma)
                   b1c  <- sum(zpsp*rs)/(n*sigma)
                   zpsi <- psi.weight(rs,ips=2,xk=1.717817)
                   a2c  <- sum(zpsi)/(n*sigma)
                   b2c  <- sum(zpsi*rs)/(n*sigma)}
if (estim=="TMLA"){sp1 <- exp(rs); sp2 <- (exp(rs)*(1+rs)-1)
                   a1c <- sum(wi*sp1)/(n*sigma)
                   b1c <- sum(wi*sp1*rs)/(n*sigma)
                   a2c <- sum(wi*sp2)/(n*sigma)
                   b2c <- sum(wi*sp2*rs)/(n*sigma)}
if (estim=="SI")  {a1c <- integrate(Pspphi.w, lower=-xk,upper=xk)$value/sigma
                   a2c <- integrate(Psiphi.w, lower=-xk,upper=xk)$value/sigma
                   b1c <- integrate(Pspzphi.w,lower=-xk,upper=xk)$value/sigma
                   b2c <- integrate(Psizphi.w,lower=-xk,upper=xk)$value/sigma}
if (estim=="TMLI"){a1c <- integrate(s1pphi.w, lower=l,  upper=u )$value/sigma
                   a2c <- integrate(s2pphi.w, lower=l,  upper=u )$value/sigma
                   b1c <- integrate(s1pzphi.w,lower=l,  upper=u )$value/sigma
                   b2c <- integrate(s2pzphi.w,lower=l,  upper=u )$value/sigma}
#cat("a1c,b1c,a2c,b2c:",round(c(a1c,b1c,a2c,b2c),4),"\n")
M  <- matrix(0,ncol=2,nrow=2)
M[1,1] <- a1c; M[1,2] <- b1c
M[2,1] <- a2c; M[2,2] <- b2c
Minv <- solve(M)
list(Minv=Minv)}

