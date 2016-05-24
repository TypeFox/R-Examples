# ex4.02.R
set.seed(123)
n.sim <- 15
aic <-matrix(, n.sim, 3)
for(i in 1:n.sim){
dri <- expression(-(x-10))
dif <- expression(2*sqrt(x)) 
sde.sim(X0=10,drift=dri, sigma=dif,N=1000,delta=0.1) -> X

b <- function(x,theta) -theta[1]*(x-theta[2])
b.x <- function(x,theta)  -theta[1]+0*x

s <- function(x,theta) theta[3]*sqrt(x)
s.x <- function(x,theta) theta[3]/(2*sqrt(x))
s.xx <- function(x,theta) -theta[3]/(4*x^1.5)
aic[i,1] <- sdeAIC(X, NULL, b, s, b.x, s.x, s.xx, guess=c(1,1,1),
            lower=rep(1e-3,3), method="L-BFGS-B")

s <- function(x,theta) sqrt(theta[3]*+theta[4]*x)
s.x <- function(x,theta) theta[4]/(2*sqrt(theta[3]+theta[4]*x))
s.xx <- function(x,theta) -theta[4]^2/(4*(theta[3]+theta[4]*x)^1.5)
aic[i,2] <- sdeAIC(X, NULL, b, s, b.x, s.x, s.xx, guess=c(1,1,1,1),
              lower=rep(1e-3,4), method="L-BFGS-B")

s <- function(x,theta) (theta[3]+theta[4]*x)^theta[5]
s.x <- function(x,theta) theta[4]*theta[5]*(theta[3]+theta[4]*x)^(-1+theta[5])
s.xx <- function(x,theta) 
     theta[4]^2*theta[5]*(theta[5]-1)*(theta[3]+theta[4]*x)^(-2+theta[5])
aic[i,3] <- sdeAIC(X, NULL, b, s, b.x, s.x, s.xx, guess=c(1,1,1,1,1),
               lower=rep(1e-3,5), method="L-BFGS-B")
}
print(aic)
table(apply(aic,1,function(x) which(x==min(x))))
