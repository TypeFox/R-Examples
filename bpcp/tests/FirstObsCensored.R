## In 1.3.0 there was an error when the first observation 
## was censored. It was a problem in the abmm function.
## fixed March 25, 2016
library(bpcp)
time<-c(1:3)
status<-c(0,1,1)
bout<-bpcp(time,status)
summary(bout)
