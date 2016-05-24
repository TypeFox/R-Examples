d <- generateFakeData()

# model fitting only
(fit <- fSAE(y0 ~ x + area2, data=d$sam, area="area"))

# model fitting and small area estimation, unit-level model
saeHB <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, silent=TRUE)
saeHB  # print a summary
EST(saeHB)  # small area estimates
SE(saeHB)  # error estimates
str(saeHB)
plot(saeHB, list(est=d$mY0), CI=2)  # compare to true population means

# unit-level model with REML model-fit instead of Bayesian approach
saeREML <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, method="REML", silent=TRUE)
plot(saeHB, saeREML)  # compare

# basic area-level model
saeA <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="area")
plot(saeHB, saeA)

# SAE estimates based on a linear unit-level model without area effects
saeL <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, method="synthetic")
plot(saeHB, saeL)

# model-based estimation of overall population mean without area effects
est.global <- fSAE(y0 ~ x + area2, data=d$sam, area=NULL, popdata=colSums(d$Xpop), method="synthetic")
EST(est.global); MSE(est.global)

# no model fitting or estimation, but return design matrix, variable of interest,
#   area indicator, area population sizes and matrix of population means
dat <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="data")
str(dat)
