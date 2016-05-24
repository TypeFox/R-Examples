## simulated data example considering Markov z's and both 
## penalized log-likelihood and Bayesian approaches. 

library(switchnpreg)

###########################################################
### setting up true parameters and generating the data  ###
###########################################################

J <- 2 ## number of states
x <- seq(from=1, to=100, by=0.5)


load(system.file(file.path('demo-support-files', 'simulated_data_fTrue.RData'),
                 package='switchnpreg'))
f_functions <- f.true

source(system.file(file.path('demo-support-files', 'generate.R'),
                 package='switchnpreg'),
       local = TRUE)
source(system.file(file.path('demo-support-files', 'initial.R'),
                 package='switchnpreg'),
       local = TRUE)

set.seed(1)
sim.data <-  make.z.y(f = f_functions,
                      sigma2 = c(0.00005,0.00005),
                      alpha = list(A=matrix(c(0.9,0.1,0.2,0.8), byrow=TRUE, nrow=2),
                                   PI=c(0.5,0.5)),
                      z.markov)

z <- sim.data$z
y <- sim.data$y

#######################################
## Penalized log-likelihood approach ##
#######################################

basis <- create.bspline.basis(range=c(min(x), max(x)), norder=4, nbasis = 40)
B <- getbasismatrix(evalarg=x, basisobj=basis, nderiv=0)
R <- getbasispenalty(basisobj=basis)            

## starting values 
initial.obj <- initial.Markov.PL.Kmeans(x,y,
                                        interval = c(log(1E-8), log(1E+6)),
                                        basis, R=R, K=3) 
f.initial <- initial.obj$f
sig2 <- initial.obj$sigma2
alpha <- list(A=matrix(c(0.5,0.5,0.5,0.5), byrow=TRUE, nrow=2),
              PI=c(0.5,0.5))

#################
## Estimation ###
#################

estimates <- switchnpreg(method = 'pl',
                         x = x, y = y,
                         f = f.initial, alpha = alpha, sig = sig2,
                         lambda = rep(6,J),
                         B=B, R=R, 
                         var.equal = TRUE,
                         z.indep = FALSE,
                         interval = c(log(1E-8), log(1E+6)),

                         eps.cv = rep(1E-6,J),
                         eps.em = c(rep(1E-4,J),rep(1E-7,J),rep(1E-3,J^2)),
                         maxit.cv=10,
                         maxit.em=100)

print('Results PL approach')
print('Transition probabilities')
print(estimates$current$alpha$A)
print("s.e.'s a_12, a_21")
print(estimates$stderr)

print(data.frame('sigma^2'=estimates$current$sigma2,
                 lambda=exp(estimates$lambda),
                 check.names=FALSE),
      digits=4)


par(mfrow=c(1,2))
plot(x, y, col=z, pch=20, main='Penalized log-likelihood approach')
matlines(x, f_functions, col=c(1,2), lty=1)
matlines(x, estimates$current$f, col=c(1,2), lty=2)
matlines(x, f.initial, col='green', lty=1)

#######################
## Bayesian approach ##
#######################

### starting values 
initial.obj <- initial.Markov.Bayes.Kmeans(x, y,
                                           Interval1=c(1,100),
                                           Interval2=c(1,100),
                                           sigma2=c(0.005,0.005),
                                           nSubInt=1, K=3)
f.initial <- initial.obj$f
sig2 <- initial.obj$sigma2
lambda <- initial.obj$lambda
alpha <-  list(A=matrix(c(0.5,0.5,0.5,0.5), byrow=TRUE, nrow=2),
               PI=c(0.5,0.5))  


#################
## Estimation ###
#################

estimates <- switchnpreg(method = 'bayes',
                         x = x, y = y,
                         f = f.initial, alpha = alpha, sig = sig2,
                         lambda = lambda,

                         var.equal = TRUE,
                         z.indep = FALSE,
                         interval = c(1, 100),
                         cov.function = switchnpreg:::covariance,

                         eps.cv = rep(1E-6,J),
                         eps.em = c(rep(1E-4,J),rep(1E-7,J),rep(1E-3,J^2)),
                         maxit.cv=10,
                         maxit.em=100)

print('Results Bayesian approach')
print('Transition probabilities')
print(estimates$current$alpha$A)
print("s.e.'s a_12, a_21")
print(estimates$stderr)

print(data.frame('sigma^2'=estimates$current$sigma2,
                 lambda=estimates$lambda,
                 check.names=FALSE),digits=4)

plot(x, y, col=z, pch=20, main='Bayesian approach')
matlines(x, f_functions, col=c(1,2), lty=1)
matlines(x, estimates$current$f, col=c(1,2), lty=2)
matlines(x, f.initial, col='green', lty=1)
