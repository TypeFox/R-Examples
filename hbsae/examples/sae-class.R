d <- generateFakeData()

# compute small area estimates
sae <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop)

coef(sae)  # fixed effects
raneff(sae)  # random effects
sv2(sae)  # between-area variance
se2(sae)  # within-area variance
cAIC(sae)  # conditional AIC
