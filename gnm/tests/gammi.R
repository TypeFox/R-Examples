library(gnm)
set.seed(1)

yield.scaled <- wheat$yield * sqrt(3/1000)
treatment <- factor(paste(wheat$tillage, wheat$summerCrop, wheat$manure,
                          wheat$N, sep = ""))
bilinear3 <- gnm(yield.scaled ~ year + treatment +
                 instances(Mult(year, treatment), 3),
                 family = gaussian, data = wheat)

print(bilinear3$deviance, digits=10)
print(bilinear3$df)
