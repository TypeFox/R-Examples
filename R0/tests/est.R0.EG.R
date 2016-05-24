#Loading package
library(R0)

## Data is taken from the paper by Nishiura for key transmission parameters of an institutional
## outbreak during 1918 influenza pandemic in Germany)

data(Germany.1918)
mGT<-generation.time("gamma", c(3, 1.5))

est.R0.EG(Germany.1918, mGT, begin=1, end=27)
## Reproduction number estimate using  Exponential Growth 
## R :  1.525895[ 1.494984 , 1.557779 ]
