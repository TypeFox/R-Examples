library("R2BayesX")

## load forest health data and tree location map
data("ForestHealth")
data("BeechBnd")

## estimate model without spatial effect
fm1 <- bayesx(defoliation ~  stand + fertilized + 
  humus + moisture + alkali + ph + soil + 
  sx(age) + sx(inclination) + sx(canopy) + 
  sx(year) + sx(elevation), family = "cumlogit", 
  method = "REML", data = ForestHealth)
summary(fm1)
plot(fm1, term = c("sx(age)", "sx(inclination)", 
  "sx(canopy)", "sx(year)", "sx(elevation)"))

## now include spatial effect
## warning: long runtime
fm2 <- bayesx(defoliation ~  stand + fertilized + 
  humus + moisture + alkali + ph + soil + 
  sx(age) + sx(inclination) + sx(canopy) + sx(year) + 
  sx(elevation) + sx(id, bs = "gk", map = BeechBnd, full = TRUE),
  family = "cumlogit", method = "REML", data = ForestHealth)
summary(fm2)
plot(fm2, term = c("sx(age)", "sx(inclination)", 
  "sx(canopy)", "sx(year)", "sx(elevation)", "sx(id)"),
  map = BeechBnd, pos = "topleft")

## compare effects for elevation and inclination
plot(c(fm1, fm2), term = c("sx(elevation)", "sx(inclination)"))
