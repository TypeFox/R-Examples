#########################################################
# Section 5.7 Monte Carlo Method for Computing Integrals
#########################################################

 p=rbeta(1000, 14.26, 23.19)
 est=mean(p^2)
 se=sd(p^2)/sqrt(1000)
 c(est,se)
