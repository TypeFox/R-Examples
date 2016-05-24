fullerT <-
function(x,p)
{
y <- as.matrix(x)
n <- nrow(y)
k <-5; dn<-0.290; an<-0; bn<-0.080; taumed <- -2.18;
x <- cbind(rep(1,n),1:n)
e <- y - x %*% solve(t(x) %*% x) %*% t(x) %*% y

ip <- as.integer( 0.5*(p+1))
sk <- (3*n-taumed^2*(ip+n))*(taumed*(k+taumed)*(ip+n))^(-1)

M1 <- fuller1(e,p)
tau <- M1$tau

if( tau > taumed) c <- -taumed + dn*(tau-taumed)
if( tau > -k & tau <= taumed) c <- ip*(tau/n) - 3/(tau+sk*(tau+k))
if( tau > -sqrt(3*n) & tau <= -k) c <- ip*(tau/n) - 3/tau
if( tau <= -sqrt(3*n)) c <- 0

b1 <- M1$coef[1,1] + c*M1$se
b2 <- min(c(b1,1))
newb <- rbind(b2,estmf(y,p,b2))
newb <- arlevel(newb,p)
return(newb)
}
