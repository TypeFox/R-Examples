#    Gaussian with identity link function, mu=beta0
#    one single series
#    Compare dthmm and mmglm1
#    R CMD BATCH --no-save dthmm-mmglm1-gaussian.R dthmm-mmglm1-gaussian.Rout.save

library(HiddenMarkov)

#------------------------------------------------------------------
#   Using dthmm

#   n = series length for each subject
#   N = number of subjects
n <- 5000
N <- 1

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

delta <- c(1, 0)

y <- dthmm(NULL, Pi=Pi, distn="norm", delta=delta, pm=list(mean=c(5, 2), sd=c(1, 1)))

y <- simulate(y, nsim=N*n, seed=10)
print(logLik(y))

tmp <- BaumWelch(y, bwcontrol(posdiff=FALSE, tol=1e-05, prt=FALSE))

print(summary(tmp))
print(logLik(tmp))


#------------------------------------------------------------------
#   Using mmglm1

glmformula <- formula(y$x ~ 1)
glmfamily <- gaussian(link="identity")
Xdesign <- model.matrix(glmformula)

beta <- matrix(c(5, 2), 
               ncol=ncol(Pi), nrow=ncol(Xdesign), byrow=TRUE)

y1 <- mmglm1(y$x, Pi, delta, glmfamily, beta, Xdesign, sigma=c(1, 1), msg=FALSE)
print(logLik(y1))

tmp1 <- BaumWelch(y1, bwcontrol(posdiff=FALSE, tol=1e-05, prt=FALSE))

print(summary(tmp1))
print(logLik(tmp1, fortran=TRUE))
print(logLik(tmp1, fortran=FALSE))

#------------------------------------------------------------------
#   Compare Models

if (abs(logLik(tmp)-logLik(tmp1)) > 1e-06)
    warning("WARNING: See tests/dthmm-mmglm1-gaussian.R, log-likelihoods are different")

if (any(Viterbi(tmp)!=Viterbi(tmp1)))
    warning("WARNING: See tests/dthmm-mmglm1-gaussian.R, Viterbi paths are different")

if (any(abs(residuals(tmp)-residuals(tmp1)) > 1e-06))
    warning("WARNING: See tests/dthmm-mmglm1-gaussian.R, residuals are different")


print(tmp$pm)
print(tmp1$beta)
print(tmp1$sigma)

print(tmp$Pi)
print(tmp1$Pi)


