#    Compare mmglm0 and mmglm1
#    Gaussian with identity link function
#    R CMD BATCH --no-save mmglm0-mmglm1-gaussian.R mmglm0-mmglm1-gaussian.Rout.save

library(HiddenMarkov)


delta <- c(0,1)

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

beta <- matrix(c(0.1, -0.1,
                 1.0,  5.0),
               byrow=TRUE, nrow=2)

sd <- c(1, 2)

#--------------------------------------------------------
#     Gaussian with identity link function
#         using mmglm0

x0 <- mmglm0(NULL, Pi, delta, family="gaussian", link="identity",
             beta=beta, sigma=sd, msg=FALSE)

x0 <- simulate(x0, nsim=5000, seed=10)

x0 <- BaumWelch(x0, bwcontrol(prt=FALSE))

print(summary(x0))

#--------------------------------------------------------
#    Now embed this data into a mmglm1 object

glmformula <- formula(y ~ x1)
glmfamily <- gaussian(link="identity")
Xdesign <- model.matrix(glmformula, data=x0$x)

x1 <- mmglm1(x0$x$y, Pi, delta, glmfamily, beta, Xdesign, sigma=sd, msg=FALSE)

x1 <- BaumWelch(x1, bwcontrol(prt=FALSE))

print(summary(x1))

#--------------------------------------------------------
#   Compare Models

if (abs(logLik(x0)-logLik(x1)) > 1e-06)
    warning("WARNING: See tests/mmglm0-mmglm1-gaussian.R, log-likelihoods are different")

if (any(Viterbi(x0)!=Viterbi(x1)))
    warning("WARNING: See tests/mmglm0-mmglm1-gaussian.R, Viterbi paths are different")

if (any(abs(residuals(x0)-residuals(x1)) > 1e-06))
    warning("WARNING: See tests/mmglm0-mmglm1-gaussian.R, residuals are different")


