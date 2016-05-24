library(quipu)

data(potato.quipu)
dat = potato.quipu

str(dat)

rquipu(dat)

rquipu(dat, layout="no text", res=c(600,400))
