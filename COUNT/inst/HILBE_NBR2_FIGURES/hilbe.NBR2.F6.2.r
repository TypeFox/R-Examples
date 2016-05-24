# hilbe.NBR2.F6.2.r
# Table 6.15 + added code
# Table of Poisson observed vs predicted mean counts and diffrence;  
#   graphic of observed vs predicted counts
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Table 6.15; Figure 6.2  
# User to amend default dataset, response, and number of counts
#
rm(list=ls())
load("c://source/medpar.RData")
mdpar <- glm(los ~ hmo+white+type2+type3, family=poisson, data=medpar)
mu <- fitted(mdpar)
avgp <- sapply(0:25, function(i) mean(exp(-mu)*(mu^i)/factorial(i)))
propObsv <- with(subset(medpar, los < 26), table(los) / nrow(medpar))
Diff <- c(0,propObsv)*100 - avgp[1:25]*100
data.frame(LOS=0:25, ObsProp=c(0,propObsv)*100, avgp*100, Diff)

plot(0:25, avgp, type="b", xlim=c(0,25),
main = "Observed vs Predicted Days",
xlab = "Days in Hospital", ylab = "Probability of LOS")
lines(0:25, c(0,propObsv), type = "b", pch = 2)
legend("topright", legend = c("Predicted Days","Observed Days"),
           lty = c(1,1), pch = c(1,2))



