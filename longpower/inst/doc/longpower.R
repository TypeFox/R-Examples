## ----knitropts, results = 'hide', echo = FALSE, message = FALSE----------
library(knitr)
opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, tidy = FALSE, comment = NA)

## ----echo=TRUE, results='hide'-------------------------------------------
library(longpower)

## ----echo=TRUE-----------------------------------------------------------
liu.liang.linear.power

## ------------------------------------------------------------------------
n = 3
t = c(0,2,5)
u = list(u1 = t, u2 = rep(0,n))
v = list(v1 = cbind(1,1,rep(0,n)),
         v2 = cbind(1,0,t))         
rho = c(0.2, 0.5, 0.8)
sigma2 = c(100, 200, 300)
tab = outer(rho, sigma2, 
      Vectorize(function(rho, sigma2){
        ceiling(diggle.linear.power(
          d=0.5,
          t=t,
          sigma2=sigma2,
          R=rho,
          alternative="one.sided",
          power=0.80)$n)}))
colnames(tab) = paste("sigma2 =", sigma2)
rownames(tab) = paste("rho =", rho)
tab

## ------------------------------------------------------------------------
# var of random intercept
sig2.i = 55
# var of random slope
sig2.s = 24
# residual var
sig2.e = 10
# covariance of slope and intercep
cov.s.i <- 0.8*sqrt(sig2.i)*sqrt(sig2.s)

cov.t <- function(t1, t2, sig2.i, sig2.s, cov.s.i){
        sig2.i + t1*t2*sig2.s + (t1+t2)*cov.s.i 
}

t = seq(0,1.5,0.25)
n = length(t)
R = outer(t, t, function(x,y){cov.t(x,y, sig2.i, sig2.s, cov.s.i)})
R = R + diag(sig2.e, n, n)
u = list(u1 = t, u2 = rep(0,n))
v = list(v1 = cbind(1,1,rep(0,n)),
         v2 = cbind(1,0,t))         

liu.liang.linear.power(d=1.5, u=u, v=v, R=R, sig.level=0.05, power=0.80)

## ------------------------------------------------------------------------
x = (rbind(1,t)%*%solve(R)%*%cbind(1,t))[2,2]
x*2*(qnorm(1-0.05/2) + qnorm(0.80))^2/1.5^2

## ------------------------------------------------------------------------
x = solve(rbind(1,t)%*%solve(R)%*%cbind(1,t))[2,2]
x*2*(qnorm(1-0.05/2) + qnorm(0.80))^2/1.5^2

## ------------------------------------------------------------------------
X = t(v[[1]])%*%solve(R)%*%v[[1]] + 
    t(v[[2]])%*%solve(R)%*%v[[2]]

Sigma1 = ((t(u[[1]])%*%solve(R)%*%t - 
           t(u[[1]])%*%solve(R)%*%v[[1]]%*%solve(X)%*%t(v[[1]])%*%solve(R)%*%t)/2)

(qnorm(1-0.05/2) + qnorm(0.80))^2/(Sigma1*(1.5)^2)

