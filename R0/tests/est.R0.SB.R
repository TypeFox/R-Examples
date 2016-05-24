#Loading package
library(R0)

## Data is taken from the paper by Nishiura for key transmission parameters of an institutional
## outbreak during 1918 influenza pandemic in Germany)

data(Germany.1918)
mGT <- generation.time("gamma", c(3,1.5))
SB <- est.R0.SB(Germany.1918, mGT)

## Results will include "most likely R(t)" (ie. the R(t) value for which the computed probability 
## is the highest), along with 95% CI, in a data.frame object
SB
# Reproduction number estimate using  Real Time Bayesian  method.
# 0 0 2.02 0.71 1.17 1.7 1.36 1.53 1.28 1.43 ...

SB$Rt.quant
# Date R.t. CI.lower. CI.upper.
# 1  1918-09-29 0.00      0.01      1.44
# 2  1918-09-30 0.00      0.01      1.42
# 3  1918-10-01 2.02      0.97      2.88
# 4  1918-10-02 0.71      0.07      1.51
# 5  1918-10-03 1.17      0.40      1.84
# 6  1918-10-04 1.70      1.09      2.24
# 7  1918-10-05 1.36      0.84      1.83
# 8  1918-10-06 1.53      1.08      1.94
# 9  1918-10-07 1.28      0.88      1.66
# 10 1918-10-08 1.43      1.08      1.77
# ...

## "Plot" will provide the most-likely R value at each time unit, along with 95CI
plot(SB)
## "Plotfit" will show the complete distribution of R for 9 time unit throughout the outbreak
plotfit(SB)
