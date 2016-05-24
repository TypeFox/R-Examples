#####################################################
# Section 5.4 A Beta-Binomial Model for Overdispersion
#####################################################

 library(LearnBayes)

 data(cancermortality)

 mycontour(betabinexch0,c(.0001,.003,1,20000),cancermortality,
   xlab="eta",ylab="K")

S=readline(prompt="Type  <Return>   to continue : ")

windows()
 mycontour(betabinexch,c(-8,-4.5,3,16.5),cancermortality,
   xlab="logit eta",ylab="log K")
