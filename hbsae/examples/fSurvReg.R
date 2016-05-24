d <- generateFakeData()

# generate design matrix, variable of interest, area indicator and population data
dat <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="data")

sae <- fSurvReg(dat$y, dat$X, dat$area, dat$Narea, dat$PopMeans)
EST(sae)  # estimates
SE(sae)  # standard errors
