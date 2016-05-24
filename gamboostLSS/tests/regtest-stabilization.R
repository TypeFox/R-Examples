###
# test functionality of stabilization

require("gamboostLSS")

## simulate Gaussian data
set.seed(0804)
x1 <- runif(1000)
x2 <- runif(1000)
x3 <- runif(1000)
x4 <- runif(1000)
mu    <- 1.5 +1 * x1 +4 * x2
sigma <- exp(1 - 0.2 * x3 -0.4 * x4)
y <- rnorm(mean = mu, sd = sigma, n = length(mu))
dat <- data.frame(x1, x2, x3, x4, y)

## do not use stabilization
m1 <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                  families=GaussianLSS(),
                  data=dat)
coef(m1)

## use stabilization via options (for backwards compatibility)
options(gamboostLSS_stab_ngrad = TRUE)
m2 <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                  families=GaussianLSS(),
                  data=dat)
coef(m2)
options(gamboostLSS_stab_ngrad = FALSE)

## now use novel interface via families:
m3 <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                  families = GaussianLSS(stabilization = "MAD"),
                  data=dat)
stopifnot(all.equal(coef(m3), coef(m2)))

## check if everything is handled correctly
GaussianLSS(stabilization = "MAD")
GaussianLSS(stabilization = "none")
res <- try(GaussianLSS(stabilization = "test"), silent = TRUE)
res


############################################################
## continue these checks for other families
dat$y <- runif(1000, min = 0.01, max = 0.99)
FAMILIES <- list(
    GaussianLSS,
    GammaLSS,
    BetaLSS,
    StudentTLSS)
for (i in 1:length(FAMILIES)) {
    m_none <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                          families = FAMILIES[[i]](stabilization = "none"),
                          data=dat)
    m_MAD <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                         families = FAMILIES[[i]](stabilization = "MAD"),
                         data=dat)
    stopifnot(tail(risk(m_none, merge = TRUE), 1) != tail(risk(m_MAD, merge = TRUE), 1))
    cat('Risks:\n  stabilization = "none":',
        tail(risk(m_none, merge = TRUE), 1),
        '\n  stabilization = "MAD":',
        tail(risk(m_MAD, merge = TRUE), 1), "\n")
}

## check as.families interface for 2:4 parametric families
dat$y <- rnorm(1000, mean = 10, sd = 1)
FAMILIES <- list(
    "NO",
    "TF")
for (i in 1:length(FAMILIES)) {
    m_none <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                          families = as.families(FAMILIES[[i]], stabilization = "none"),
                          data=dat)
    m_MAD <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                         families = as.families(FAMILIES[[i]], stabilization = "MAD"),
                         data=dat)
    cat('Risks:\n  stabilization = "none":',
        tail(risk(m_none, merge = TRUE), 1),
        '\n  stabilization = "MAD":',
        tail(risk(m_MAD, merge = TRUE), 1), "\n")
}

FAMILIES <- list("BCT")
require("gamlss.dist")
dat$y <- rBCT(1000, mu = 100, sigma = 0.1, nu = 0, tau = 2)
for (i in 1:length(FAMILIES)) {
    m_none <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                          families = as.families(FAMILIES[[i]], stabilization = "none"),
                          data=dat)
    m_MAD <- try(glmboostLSS(y ~ x1 + x2 + x3 + x4,
                             families = as.families(FAMILIES[[i]], stabilization = "MAD"),
                             data=dat), silent = TRUE)
    if (inherits(m_MAD, "try-error")) {
        warning("BCT cannot be fitted with stabilization", immediate. = TRUE)
        break
    }
    cat('Risks:\n  stabilization = "none":',
        tail(risk(m_none, merge = TRUE), 1),
        '\n  stabilization = "MAD":',
        tail(risk(m_MAD, merge = TRUE), 1), "\n")
}

## check survival families
dat$zens <- sample(c(0, 1), 1000, replace = TRUE)
FAMILIES <- list(
    LogNormalLSS,
    WeibullLSS,
    LogLogLSS)
families <- c("LogNormalLSS", "WeibullLSS", "LogLogLSS")
require(survival)
for (i in 1:length(FAMILIES)) {
    cat(families[i], "\n\n")
    m_none <- glmboostLSS(Surv(y, zens) ~ x1 + x2 + x3 + x4,
                          families = FAMILIES[[i]](stabilization = "none"),
                          data=dat)
    m_MAD <- try(glmboostLSS(Surv(y, zens) ~ x1 + x2 + x3 + x4,
                             families = FAMILIES[[i]](stabilization = "MAD"),
                             data=dat), silent = TRUE)

    if (inherits(m_MAD, "try-error")) {
        warning(families[i], "cannot be fitted with stabilization", immediate. = TRUE)
        break
    }
    cat('Risks:\n  stabilization = "none":',
        tail(risk(m_none, merge = TRUE), 1),
        '\n  stabilization = "MAD":',
        tail(risk(m_MAD, merge = TRUE), 1), "\n")
}

## check count data
dat$y <- rnbinom(1000, size=10, mu=5)
FAMILIES <- list(
    NBinomialLSS,
    ZIPoLSS,
    ZINBLSS)
for (i in 1:length(FAMILIES)) {
    m_none <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                          families = FAMILIES[[i]](stabilization = "none"),
                          data=dat)
    m_MAD <- glmboostLSS(y ~ x1 + x2 + x3 + x4,
                         families = FAMILIES[[i]](stabilization = "MAD"),
                         data=dat)
    cat('Risks:\n  stabilization = "none":',
        tail(risk(m_none, merge = TRUE), 1),
        '\n  stabilization = "MAD":',
        tail(risk(m_MAD, merge = TRUE), 1), "\n")
}
