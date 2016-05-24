d <- generateFakeData()

# generate design matrix, variable of interest, area indicator and population data
dat <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="data")

# compute small area estimates based on the basic unit-level model
sae <- fSAE.Unit(dat$y, dat$X, dat$area, dat$Narea, dat$PopMeans)
EST(sae)  # estimates
SE(sae)  # standard errors
