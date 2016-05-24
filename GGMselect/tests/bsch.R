library(GGMselect)
attach(system.file("extdata", "G3simone.Rdata", package="GGMselect"))
attach(system.file("extdata", "genes.lm.Rdata", package="GGMselect"))

print(selectMyFam(genes.lm, MyFamily = G3simone, K=10))
