# hilbe.NBR2.F6.2alt.r
# Table 6.15 + added code
# Table of Poisson observed vs predicted mean counts and diffrence;  
#   graphic of observed vs predicted counts
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Table 6.15; Figure 6.2  Alternative
# User to amend default dataset, response, and number of counts
#
load("c://source/medpar.RData") 
mdpar <- glm(los ~ hmo+white+type2+type3,family=poisson, data=medpar) 
mu <- fitted.values(mdpar) 
p <- NULL 
avgp <- NULL
for (i in 0:25) { 
 p[[i+1]] <- exp(-mu)*(mu^i)/factorial(i) 
 avgp[i+1] <- mean(p[[i+1]])
} 
nCases <- dim(medpar) 
n<- NULL 
propObs<- NULL 
probFit<- NULL 
yFitMean<- NULL 
for (i in 0:25) {                  #possible values for LOS 
bLos<- medpar$los==i               #selector for los=i 
n[i+1]<- sum(bLos) #number for los=i 
propObs[i+1]<- n[i+1]/nCases[1]    #observed proportion for LOS=i 
} 
Diff <- propObs*100 - avgp*100
data.frame(LOS=0:25, ObsProp=propObs*100, avgp*100, Diff)


