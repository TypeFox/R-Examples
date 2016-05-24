# hilbe.NBR2.F9.2alt.r
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Figure 9.2  Table Creation - no graphic  Alternative  
# 
rm(list=ls())
load("c://source/affairs.RData") 
affmodel <- glm.nb(naffairs~ kids + avgmarr + hapavg +  vryhap + notrel + slghtrel + smerel + vryrel + yrsmarr3 + yrsmarr4 + yrsmarr5 + yrsmarr6, data=affairs)
mu <- fitted.values(affmodel) 
avgp <- sapply(0:13, function(i) mean(exp(-mu)*(mu^i)/factorial(i)))
propObsv <- with(subset(affairs, naffairs < 14), table(naffairs) / nrow(affairs))
Diff <- c(0,propObsv)*100 - avgp[1:25]*100
data.frame(LOS=0:13, ObsProp=c(0,propObsv)*100, avgp*100, Diff)