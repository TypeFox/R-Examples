
library(RApiSerialize)

data(trees)
fit <- lm(log(Girth) ~ log(Volume) + log(Height), trees)

## serialize and use R's unserialize
identical(unserialize(serializeToRaw(fit)), fit)
## serialize and use our unserialize
identical(unserializeFromRaw(serializeToRaw(fit)), fit)
## R's serialize and our unserialize
identical(unserializeFromRaw(serialize(fit, NULL)), fit)
## R's serialize and R's unserialize (doh)
identical(unserialize(serialize(fit, NULL)), fit)
