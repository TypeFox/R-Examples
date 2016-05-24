# example 1: GasolineYield
library(zoib)
data("GasolineYield", package = "zoib")

# zoib: fixed
eg.fixed <- zoib(yield ~ temp + as.factor(batch)| 1, data=GasolineYield,
                 joint = FALSE,  random = 0, EUID = 1:nrow(d),
                 zero.inflation = FALSE, one.inflation = FALSE,
                 n.iter = 50, n.thin = 2, n.burn=1)
sample1 <- eg.fixed$coeff


# zoib: random
eg.random <- zoib(yield ~ temp | 1 | 1, data=GasolineYield,
                  joint = FALSE, random=1, EUID=GasolineYield$batch,
                  zero.inflation = FALSE, one.inflation = FALSE,
                  n.iter=50, n.thin=2, n.burn=1)
sample2 <- eg.random$coeff


