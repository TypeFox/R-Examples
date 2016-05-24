###################################################
### chunk number 1: setup
###################################################
options(prompt = "R> ", continue = "+  ", width = 64,
  digits = 4, show.signif.stars = FALSE, useFancyQuotes = FALSE)

options(SweaveHooks = list(onefig =   function() {par(mfrow = c(1,1))},
                           twofig =   function() {par(mfrow = c(1,2))},                           
                           threefig = function() {par(mfrow = c(1,3))},
                           fourfig =  function() {par(mfrow = c(2,2))},
			   sixfig =   function() {par(mfrow = c(3,2))}))

library("AER")

set.seed(1071)


###################################################
### chunk number 2: DGP
###################################################
dgp <- function(nobs = 15, model = c("trend", "dynamic"),
  corr = 0, coef = c(0.25, -0.75), sd = 1)
{
  model <- match.arg(model)
  coef <- rep(coef, length.out = 2)

  err <- as.vector(filter(rnorm(nobs, sd = sd), corr,
    method = "recursive"))
  if(model == "trend") {
    x <- 1:nobs
    y <- coef[1] + coef[2] * x + err
  } else {
    y <- rep(NA, nobs)
    y[1] <- coef[1] + err[1]
    for(i in 2:nobs)
      y[i] <- coef[1] + coef[2] * y[i-1] + err[i]
    x <- c(0, y[1:(nobs-1)])
  }
  return(data.frame(y = y, x = x))
}


###################################################
### chunk number 3: simpower
###################################################
simpower <- function(nrep = 100, size = 0.05, ...)
{
  pval <- matrix(rep(NA, 2 * nrep), ncol = 2)
  colnames(pval) <- c("dwtest", "bgtest")
  for(i in 1:nrep) {
    dat <- dgp(...)
    pval[i,1] <- dwtest(y ~ x, data = dat,
      alternative = "two.sided")$p.value
    pval[i,2] <- bgtest(y ~ x, data = dat)$p.value
  }
  return(colMeans(pval < size))
}


###################################################
### chunk number 4: simulation-function
###################################################
simulation <- function(corr = c(0, 0.2, 0.4, 0.6, 0.8,
  0.9, 0.95, 0.99), nobs = c(15, 30, 50),
  model = c("trend", "dynamic"), ...)
{
  prs <- expand.grid(corr = corr, nobs = nobs, model = model)
  nprs <- nrow(prs)

  pow <- matrix(rep(NA, 2 * nprs), ncol = 2)
  for(i in 1:nprs) pow[i,] <- simpower(corr = prs[i,1],
    nobs = prs[i,2], model = as.character(prs[i,3]), ...)

  rval <- rbind(prs, prs)
  rval$test <- factor(rep(1:2, c(nprs, nprs)),
    labels = c("dwtest", "bgtest"))
  rval$power <- c(pow[,1], pow[,2])
  rval$nobs <- factor(rval$nobs)
  return(rval)
}


###################################################
### chunk number 5: simulation
###################################################
set.seed(123)
psim <- simulation()


###################################################
### chunk number 6: simulation-table
###################################################
tab <- xtabs(power ~ corr + test + model + nobs, data = psim)
ftable(tab, row.vars = c("model", "nobs", "test"),
  col.vars = "corr")


###################################################
### chunk number 7: simulation-visualization
###################################################
library("lattice")
xyplot(power ~ corr | model + nobs, groups = ~ test,
  data = psim, type = "b")


###################################################
### chunk number 8: simulation-visualization1
###################################################
library("lattice")
trellis.par.set(theme = canonical.theme(color = FALSE))
print(xyplot(power ~ corr | model + nobs, groups = ~ test, data = psim, type = "b"))


###################################################
### chunk number 9: journals-lm
###################################################
data("Journals")
journals <- Journals[, c("subs", "price")]
journals$citeprice <- Journals$price/Journals$citations
jour_lm <- lm(log(subs) ~ log(citeprice), data = journals)


###################################################
### chunk number 10: journals-residuals-based-resampling-unused eval=FALSE
###################################################
## refit <- function(data, i) {
##   d <- data
##   d$subs <- exp(d$fitted + d$res[i])
##   coef(lm(log(subs) ~ log(citeprice), data = d))
## }


###################################################
### chunk number 11: journals-case-based-resampling
###################################################
refit <- function(data, i)
  coef(lm(log(subs) ~ log(citeprice), data = data[i,]))


###################################################
### chunk number 12: journals-boot
###################################################
library("boot")
set.seed(123)
jour_boot <- boot(journals, refit, R = 999)


###################################################
### chunk number 13: journals-boot-print
###################################################
jour_boot


###################################################
### chunk number 14: journals-lm-coeftest
###################################################
coeftest(jour_lm)


###################################################
### chunk number 15: journals-boot-ci
###################################################
boot.ci(jour_boot, index = 2, type = "basic")


###################################################
### chunk number 16: journals-lm-ci
###################################################
confint(jour_lm,  parm = 2)


###################################################
### chunk number 17: ml-loglik
###################################################
data("Equipment", package = "AER")

nlogL <- function(par) {
  beta <- par[1:3]
  theta <- par[4]
  sigma2 <- par[5]

  Y <- with(Equipment, valueadded/firms)
  K <- with(Equipment, capital/firms)
  L <- with(Equipment, labor/firms)

  rhs <- beta[1] + beta[2] * log(K) + beta[3] * log(L)
  lhs <- log(Y) + theta * Y

  rval <- sum(log(1 + theta * Y) - log(Y) +
    dnorm(lhs, mean = rhs, sd = sqrt(sigma2), log = TRUE))
  return(-rval)
}


###################################################
### chunk number 18: ml-0
###################################################
fm0 <- lm(log(valueadded/firms) ~ log(capital/firms) +
  log(labor/firms), data = Equipment)


###################################################
### chunk number 19: ml-0-coef
###################################################
par0 <- as.vector(c(coef(fm0), 0, mean(residuals(fm0)^2)))


###################################################
### chunk number 20: ml-optim
###################################################
opt <- optim(par0, nlogL, hessian = TRUE)


###################################################
### chunk number 21: ml-optim-output
###################################################
opt$par
sqrt(diag(solve(opt$hessian)))[1:4]
-opt$value


###################################################
### chunk number 22: Sweave eval=FALSE
###################################################
## Sweave("Sweave-journals.Rnw")


###################################################
### chunk number 23: Stangle eval=FALSE
###################################################
## Stangle("Sweave-journals.Rnw")


###################################################
### chunk number 24: texi2dvi eval=FALSE
###################################################
## texi2dvi("Sweave-journals.tex", pdf = TRUE)


###################################################
### chunk number 25: vignette eval=FALSE
###################################################
## vignette("Sweave-journals", package = "AER")


