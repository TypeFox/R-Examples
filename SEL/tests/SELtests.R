library(SEL)
### Check calculation of the entropy term
foo <- function(x, fit){
  predict(fit, newdata = x)^2
}
### Case 1
x <- c(2,5,7,9)
y <- c(0.05, 0.2, 0.8, 0.9)
fit1 <- SEL(x, y, gamma = 1, bounds = c(0,10))
as.numeric(fit1$coefs%*%fit1$Omega%*%fit1$coefs)
integrate(foo, 0, 10, fit=fit1)

### Case 2
x <- c(-0.9, -0.5, -0.1)
y <- c(0.1, 0.8, 0.9)
fit2 <- SEL(x, y, Delta = 0.03, bounds = c(-1,0), N = 10)
as.numeric(fit2$coefs%*%fit2$Omega%*%fit2$coefs)
integrate(foo, -1, 0, fit=fit2)

### Case 3
x <- seq(1e5, 2e5, length = 10)
y <- c(seq(0.1, 0.4, length = 5), seq(0.65,0.9,length=5))
fit3 <- SEL(x, y, bounds = c(5e4, 2.5e5), N=0, d=20, gamma = 2,
            constr = "unimodal", mode = 1.5e5)
as.numeric(fit3$coefs%*%fit3$Omega%*%fit3$coefs)
integrate(foo, 5e4, 2.5e5, fit=fit3)

### Case 4
pistar <- function(x) dbeta(x, 1.5, 1.5)
x <- c(0.1, 0.5, 0.9)
y <- c(0.1, 0.8, 0.9)
fit2 <- SEL(x, y, Delta = 0.03, bounds = c(0,1), N = 10, pistar=pistar)
fit2
sq <- seq(0, 1, length=2001)
mean((predict(fit2, newdata=sq)-dbeta(sq, 1.5, 1.5))^2)
