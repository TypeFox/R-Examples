#################################################
# Section 8.4  A Two-Sided Test of a Normal Mean
#################################################

library(LearnBayes)

 weights=c(182, 172, 173, 176, 176, 180, 173, 174, 179, 175)
 data=c(mean(weights),length(weights),3)
 t=c(.5,1,2,4,8)
 mnormt.twosided(170,.5,t,data)
