#Loading package
library(R0)

## Data is taken from paper by Nishiura for key transmission parameters of an institutional
## outbreak during the 1918 influenza pandemic in Germany)

data(Germany.1918)
mGT<-generation.time("gamma", c(2.45, 1.38))
est.R0.ML(Germany.1918, mGT, begin=1, end=27, range=c(0.01,50))
# Reproduction number estimate using  Maximum Likelihood  method.
# R :  1.307222[ 1.236913 , 1.380156 ]

res=est.R0.ML(Germany.1918, mGT, begin=1, end=27, range=c(0.01,50))
plot(res)

## no change in R with varying range 
## (dates here are the same index as before. Just to illustrate different use)
est.R0.ML(Germany.1918, mGT, begin="1918-09-29", end="1918-10-25", range=c(0.01,100))
# Reproduction number estimate using  Maximum Likelihood  method.
# R :  1.307249[ 1.236913 , 1.380185 ]


