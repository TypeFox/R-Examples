d <- generateFakeData()

# compute small area estimates based on area-level and unit-level models
saeArea <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="area",
  silent=TRUE, keep.data=TRUE)
saeUnit <- fSAE(y0 ~ x + area2, data=d$sam, area="area", popdata=d$Xpop, type="unit",
  silent=TRUE, keep.data=TRUE)

# compare area and unit-level models based on area-level cross-validation criterion
CVarea(saeArea, K=10, seed=1)  # 10-fold CV for area-level model
CVarea(saeUnit, K=10, seed=1)  # 10-fold CV for unit-level model
