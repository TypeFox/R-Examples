##########################################################
###  Motorcycle data example considering that the data ###
###  came from J=3 simulated accidents                 ###
##########################################################

library(switchnpreg)

initial.shift <-  c(30, 0, -30)  ## for initial functions, must be of length J
J <- length(initial.shift)

x <- MASS::mcycle$times
set.seed(30)
x[duplicated(x)] <- round(jitter(x[duplicated(x)]),3)

y <- MASS::mcycle$accel

n <- length(y)
    
spline_fit <- smooth.spline(x, y)

#######################################
## setting up the initial functions  ##
#######################################
    
f.initial <- matrix(spline_fit$y, nrow = n, ncol = J)
f.initial <- t(apply(f.initial, 1, `+`, initial.shift))

    
############################################################
## estimation using the penalized log-likelihood approach ##
############################################################
                          
basis <- create.bspline.basis(range=c(min(x), max(x)), norder=4, nbasis = 40)
B <- getbasismatrix(evalarg=x, basisobj=basis, nderiv=0)
R <- getbasispenalty(basisobj=basis)


## initial values PL approach
lambda <- rep(.5, J) ## in log scale
sig2 <- rep((sum((y-predict(spline_fit, x)$y)^2) / (n - spline_fit$df))/J, J)
alpha <- rep(1, J) / J

estimates <- switchnpreg(method = 'pl',
                         x = x, y = y,
                         f = f.initial, alpha = alpha, sig = sig2,
                         lambda = rep(.5,J),
                         B=B, R=R,
                         var.equal = FALSE,
                         interval = c(log(10^(-4)), log(10^3)),

                         eps.cv = rep(1E-1,J),
                         eps.em = rep(c(1E-1, 1E-2, 1E-3), each = J),
                         maxit.cv=10,
                         maxit.em=100)

print('Results PL approach')
print(round(data.frame('sigma^2'=estimates$current$sigma2,
                       p=estimates$current$alpha,
                       'stderr p'= estimates$stderr,
                       lambda=exp(estimates$lambda),
                       check.names=FALSE),
            digits=3))

par(mfrow=c(1,2))
plot(x, y, ylim=c(-150,90),
     ylab='Head acceleration',
     xlab='Time',
     main='Penalized log-likelihood approach, J=3')
matlines(x, estimates$current$f, type='l', lty=1, col=1:J)
matlines(sort(x), f.initial, lty=2, col='gray')


############################################
## estimation using the Bayesian approach ##
############################################

## initial values Bayesian approach
lambda <- rep(5, J)
sig2 <- rep(10, J)
alpha <- rep(1, J) / J
u <- var(y) - sum((y-predict(spline_fit, x)$y)^2) / (n - spline_fit$df)


estimates <- switchnpreg(method = 'bayes',
                         x = x, y = y,
                         f = f.initial, alpha = alpha, sig = sig2,
                         lambda = lambda,

                         var.equal = FALSE,
                         interval = c(1, 100),
                         cov.function = covariance_fixed_u,
                         u = u,

                         eps.cv = rep(1E-1,J),
                         eps.em = rep(c(1E-1, 1E-2, 1E-3), each = J),
                         maxit.cv=10,
                         maxit.em=100)


print('Results Bayesian approach')
print(round(data.frame('sigma^2'=estimates$current$sigma2,
                       p=estimates$current$alpha,
                       'stderr p'= estimates$stderr,
                       lambda=estimates$lambda,
                       check.names=FALSE),
            digits=3))


plot(x, y,
     ylim=c(-150,90),
     ylab='Head acceleration',
     xlab='Time',
     main='Bayesian approach, J=3')
matlines(x, estimates$current$f, type='l', lty=1, col=1:J)
matlines(sort(x), f.initial, lty=2, col='gray')
