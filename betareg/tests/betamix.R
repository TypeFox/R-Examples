## package and data
library("betareg")
data("ReadingSkills", package = "betareg")

## beta regression mixture model
set.seed(1071)
rs_mix <- betamix(accuracy ~ iq, data = ReadingSkills, k = 3,
  extra_components = extraComponent(type = "uniform",
    coef = 0.99, delta = 0.01), nstart = 10)

## fitted model
print(rs_mix)
summary(rs_mix)

## further methods
table(clusters(rs_mix), ReadingSkills$dyslexia)
posterior(rs_mix)
