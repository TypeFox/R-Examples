#    Compare mmglm0 and mmglm1
#    Gaussian with identity link function
#    R CMD BATCH --no-save mmglm0-mmglm1-binomial.R mmglm0-mmglm1-binomial.Rout.save

library(HiddenMarkov)


delta <- c(0,1)

Pi <- matrix(c(0.8, 0.2,
               0.3, 0.7),
             byrow=TRUE, nrow=2)

beta <- matrix(c(0.1, -0.1,
                 1.0,  5.0),
               byrow=TRUE, nrow=2)

sd <- c(1, 2)

n <- 5000

#   Use different numbers of Bernoulli trials
set.seed(5)
x <- list(size=rpois(n, 10)+1)

#--------------------------------------------------------
#     Gaussian with identity link function
#         using mmglm0

x0 <- mmglm0(x, Pi, delta, family="binomial", link="logit",
             beta=beta, sigma=sd, msg=FALSE)

x0 <- simulate(x0, nsim=n, seed=10)

x0 <- BaumWelch(x0, bwcontrol(prt=FALSE))

print(summary(x0))

#--------------------------------------------------------
#    Now embed this data into a mmglm1 object

glmformula <- formula(y ~ x1)
glmfamily <- binomial(link="logit")
Xdesign <- model.matrix(glmformula, data=x0$x)

x1 <- mmglm1(x0$x$y, Pi, delta, glmfamily, beta, Xdesign, sigma=sd,
             size=x$size, msg=FALSE)

x1 <- BaumWelch(x1, bwcontrol(prt=FALSE))

print(summary(x1))

#--------------------------------------------------------
#   Compare Models

if (abs(logLik(x0)-logLik(x1)) > 1e-06)
    warning("WARNING: See tests/mmglm0-mmglm1.R-binomial, log-likelihoods are different")

if (any(Viterbi(x0)!=Viterbi(x1)))
    warning("WARNING: See tests/mmglm0-mmglm1-binomial.R, Viterbi paths are different")

if (any(abs(residuals(x0)-residuals(x1)) > 1e-06))
    warning("WARNING: See tests/mmglm0-mmglm1-binomial.R, residuals are different")


