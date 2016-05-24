load("../testdata.rda")

###############################################################################
## Tests compare ordinalNet results with other packages in special cases where
## comparisons are possible:
## e.g. unpenalized ordinal logistic regression can be compared with MASS::polr
## e.g. penalized binary logistic regression can be compared with glmnet and penalized
##
## Data and results from other packages are stored in testdata.rda so that the
## tests can be run without package and version dependencies.
###############################################################################

###############################################################################
## Make multinomial matrices and link function
###############################################################################
makeMultinomialXLS <- function(x, nLev)
{
    splitX <- split(x, row(x))
    xLS <- lapply(splitX, function(x) {
        xx <- c(1, x)
        m1 <- lapply(1:(k-1), function(i) do.call(rbind, c(rep(list(0), i-1), list(xx), rep(list(0), k-1-i))))
        m2 <- list(do.call(rbind, rep(list(-xx), k-1)))
        mat <- do.call(cbind, c(m1, m2))
        mat
    })
    xLS
}
xLS <- makeMultinomialXLS(xScaled, k)

yMat <- matrix(0, nrow=n, ncol=k)
yMat[cbind(1:n, y)] <- 1

makeMultinomialLink <- function()
{
    lf <- list()
    lf$g <- function(p) sapply(p, function(pp) log(pp / (1-sum(p))))
    lf$h <- function(eta) sapply(eta, function(ee) exp(ee) / (1 + sum(exp(eta))))
    lf$getQ <- function(eta) {
        p <- lf$h(eta)
        q <- -tcrossprod(p)
        diag(q) <- p * (1-p)
        q
    }
    lf
}

multinomialLink <- makeMultinomialLink()

###############################################################################

test_that("Unpenalized ordinal logit matches MASS::polr", {
#     p2 <- polr(y~x)
#     p2_coef <- c(p2$zeta, -p2$coef)
    o2 <- ordinalNet(x, y, lambdaVals=0)
    expect_equal(coef(o2), p2_coef, check.attributes=F, tolerance=1e-3)
})

test_that("Elastic net binary logistic regression matches glmnet", {
#     g3 <- glmnet(x, yy, family="binomial", alpha=.5, lambda=.1)
#     g3_coef <- c(-g3$a0, -as.vector(g3$beta))
    o3 <- ordinalNet(x, yy, alpha=.5, lambdaVals=.1, epsIn=1e-8, epsOut=1e-8)
    expect_equal(coef(o3), g3_coef, check.attributes=F, tolerance=1e-3)
})

test_that("Elastic net binary logistic regression matches penalized", {
#     p3 <- penalized(yy, x, model="logistic", lambda1=.5*.1*n, lambda2=.5*.1*n, standardize=T, trace=F)
#     p3_coef <- c(-p3@unpenalized, -p3@penalized)
    o3 <- ordinalNet(x, yy, alpha=.5, lambdaVals=.1, epsIn=1e-8, epsOut=1e-8)
    expect_equal(coef(o3), p3_coef, check.attributes=F, tolerance=1e-3)
})

test_that("Elastic net binary logistic regression with positve constraints and some unpenalized terms matches penalized", {
#     p4 <- penalized(yy, penalized=-x[,-1], unpenalized=cbind(Intercept=-1, -x[,1]), model="logistic",
#                     lambda1=.5*.1*n, lambda2=.5*.1*n, positive=T, standardize=F, trace=F)
#     p4_coef <-  c(p4@unpenalized, p4@penalized)
    o4 <- ordinalNet(x, yy, standardize=F, lambdaVals=.1, alpha=.5,
                     penalizeID=c(F, rep(T, p-1)), positiveID=rep(T, p), epsIn=1e-8, epsOut=1e-8)
    expect_equal(coef(o4), p4_coef, check.attributes=F, tolerance=1e-3)
})

test_that("Best AIC ordinal logit matches ordinalgmifs", {
#     og6 <- ordinal.gmifs(y ~ 1, x=names(data.frame(x)), data.frame(x), eps=.01, scale=F)
#     og6_coef <- coef(og6)
    o6 <- ordinalNet(x, y, standardize=F)
    expect_equal(coef(o6), og6_coef, check.attributes=F, tolerance=.05)
})

test_that("Elastic net multinomial logistic regression matches glmnet", {
#     g7 <- glmnet(xScaled, y, alpha=.5, lambda=.1, family="multinomial", standardize=F)
#     g7_coef <- as.vector(Reduce(rbind, coef(g7)))
    m7 <- mirlsNet(xLS, yMat, alpha=.5, lambdaVals=.1, linkfun=multinomialLink,
                   penalizeID=rep(c(F, rep(T, p)), k), betaStart=rep(0, (p+1)*k),
                   epsIn=1e-8, epsOut=1e-8)
    shift <- coef(m7)[1] - g7_coef[1]
    g7_coef[seq(1, (p+1)*(k-1)+1, length.out=k)] <- g7_coef[seq(1, (p+1)*(k-1)+1, length.out=k)] + shift
    expect_equal(coef(m7), g7_coef, check.attributes=F, tolerance=1e-3)
})

test_that("Unpenalized multinomial logistic regression matches nnet::multinom", {
#     mu8 <- multinom(y ~ xScaled)
#     mu8_coef <- c(t(coef(mu8)))
    m8 <- mirlsNet(xLS, yMat, lambdaVals=0, linkfun=multinomialLink,
                   penalizeID=rep(c(F, rep(T, p)), k), betaStart=rep(0, (p+1)*k),
                   epsIn=1e-8, epsOut=1e-8)
    m8_coef <- coef(m8)[-(1:(p+1))] - rep(coef(m8)[1:(p+1)], k-1)
    expect_equal(m8_coef, mu8_coef, check.attributes=F, tolerance=1e-3)
})
