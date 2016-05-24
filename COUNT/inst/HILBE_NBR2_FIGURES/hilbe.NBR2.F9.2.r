# hilbe.NBR2.F9.2.r
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Figure 9.2     user may amend as required for own model
# 
load("c://source/affairs.RData") 
affmodelnb <- glm.nb(naffairs~ kids + avgmarr + hapavg +  vryhap + notrel + slghtrel + smerel + vryrel + yrsmarr3 + yrsmarr4 + yrsmarr5 + yrsmarr6, data=affairs) 
mu <- fitted.values(affmodelnb) 
p <- NULL 
avgp <- NULL
for (i in 0:15) { 
 p[[i+1]] <- exp(-mu)*(mu^i)/factorial(i) 
 avgp[i+1] <- mean(p[[i+1]])
} 
nCases <- dim(affairs) 
n<- NULL 
propObs<- NULL 
probFit<- NULL 
yFitMean<- NULL 
for (i in 0:15) {                        #possible values for naffairs
bLos<- affairs$naffairs==i               #selector for naffairs=i 
n[i+1]<- sum(bLos)                       #number for naffairs=i 
propObs[i+1]<- n[i+1]/nCases[1]          #observed proportion for naffairs=i 
} 
Diff <- propObs*100 - avgp*100
data.frame(LOS=0:15, ObsProp=propObs*100, avgp*100, Diff)

plot(0:15, avgp, type="b", xlim=c(0,15), 
      main = "Observed vs Predicted Affairs", 
      xlab = "Number Affairs", ylab = "Probability of naffairs")
lines(0:15, propObs, type = "b", pch = 2)
legend("topright", legend = c("Predicted Affairs","Observed Affairs"),
        lty = c(1,1), pch = c(1,2))  














