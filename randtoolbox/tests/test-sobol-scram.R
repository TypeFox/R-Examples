library(randtoolbox)

#5742


umat<- sobol(n=2^15,dim=12,scrambling=3,seed=1776)
umat[10253:10258,1]
warnings()
