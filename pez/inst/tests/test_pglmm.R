require(testthat)
require(pez)

test_that("pglmm",{
    nspp <- 15
    nsite <- 10
    sd.resid <- 0
    beta0 <- 0
    beta1 <- 0
    sd.B0 <- 1
    sd.B1 <- 1
    phy <- rtree(n = nspp)
    phy <- compute.brlen(phy, method = "Grafen", power = 0.5)
    Vphy <- vcv(phy)
    Vphy <- Vphy/(det(Vphy)^(1/nspp))
    X <- matrix(1:nsite, nrow = 1, ncol = nsite)
    X <- (X - mean(X))/sd(X)
    iD <- t(chol(Vphy))
    b0 <- beta0 + iD %*% rnorm(nspp, sd = sd.B0)
    b1 <- beta1 + iD %*% rnorm(nspp, sd = sd.B1)
    y <- matrix(outer(b0, array(1, dim = c(1, nsite))), nrow = nspp,
                ncol = nsite) + matrix(outer(b1, X), nrow = nspp, ncol = nsite)
    e <- rnorm(nspp * nsite, sd = sd.resid)
    y <- y + matrix(e, nrow = nspp, ncol = nsite)
    y <- matrix(y, nrow = nspp * nsite, ncol = 1)
    Y <- rbinom(n = length(y), size = 1, prob = exp(y)/(1 + exp(y)))
    Y <- matrix(Y, nrow = nspp, ncol = nsite)
    rownames(Y) <- 1:nspp
    colnames(Y) <- 1:nsite
    YY <- matrix(Y, nrow = nspp * nsite, ncol = 1)
    XX <- matrix(kronecker(X, matrix(1, nrow = nspp, ncol = 1)), nrow =
                     nspp * nsite, ncol = 1)
    site <- matrix(kronecker(1:nsite, matrix(1, nrow = nspp, ncol =
                                                 1)), nrow = nspp * nsite, ncol = 1)
    sp <- matrix(kronecker(matrix(1, nrow = nsite, ncol = 1), 1:nspp),
                 nrow = nspp * nsite, ncol = 1)
    dat <- data.frame(Y = YY, X = XX, site = as.factor(site), sp = as.factor(sp))
    re.1 <- list(1, sp = dat$sp, covar = diag(nspp))
    re.2 <- list(1, sp = dat$sp, covar = Vphy)
    re.3 <- list(dat$X, sp = dat$sp, covar = diag(nspp))
    re.4 <- list(dat$X, sp = dat$sp, covar = Vphy)
    re.site <- list(1, site = dat$site, covar = diag(nsite))
    model <- communityPGLMM(Y~X,data=dat,family="binomial",sp=dat$sp,site=dat$site,random.effects=list(re.1,re.2,re.3,re.4),REML=TRUE,verbose=FALSE)
    expect_that(names(model),equals(c("formula","data","family","random.effects","B","B.se","B.cov","B.zscore","B.pvalue","ss","s2n","s2r","s2resid","logLik","AIC","BIC","REML","s2.init","B.init","Y","X","H","iV","mu","nested","sp","site","Zt","St","convcode","niter")))
    norm <- communityPGLMM(Y~X,data=dat,sp=dat$sp,site=dat$site,random.effects=list(re.1,re.2,re.3,re.4),REML=TRUE,verbose=FALSE)
    expect_that(names(model), equals(names(norm)))
})
