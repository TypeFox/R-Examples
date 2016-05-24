library(gnm)
set.seed(1)

biplotModel <- gnm(y ~ -1 + instances(Mult(site, variety), 2),
                   family = wedderburn, data = barley)

print(biplotModel$deviance, digits = 10)
print(biplotModel$df)

