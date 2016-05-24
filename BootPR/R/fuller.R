fuller <-
function(x,p)
{
y <- as.matrix(x)
n <- nrow(y)
x <- matrix(1,nrow=n)
e <- y - x %*% solve(t(x) %*% x) %*% t(x) %*% y

M1 <- fuller1(e,p)
tau <- M1$tau

r <- 1; K <- 5; dn <- 0.1111; an <- 0.0467; bn <- 0.0477; taumed <- -1.57
k1 <- (r+1)*n/(0.5*(p+1))
k2 <- ( (1+(0.5*(p+1)/n)) * taumed*(taumed-K))^(-1)*( (r+1)-(0.5*(p+1)/n)*taumed^2) 

if(tau <= -sqrt(k1)) c <- 0
if(tau > -sqrt(k1) & tau <= -K) c <- 0.5*(p+1)/n*tau - (r+1)/tau
if(tau > -K & tau <= taumed) c <- 0.5*(p+1)/n*tau - (r+1)*(tau+k2*(tau+K))^(-1)
if(tau >= taumed) c <- -taumed + dn*(tau-taumed)

b1 <- M1$coef[1,1] + c*M1$se
b2 <- min(c(b1,1))
{
if( p ==1 & b2 == 1) newb <- rbind(b2,0)
else
{
newb <- rbind(b2,estmf0(y,p,b2))
newb <- arlevel(newb,p)}
}
return(newb)
}
