# hilbe.NBR2.F9.5.r
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# observed vs predicted number of afffairs
# Figure 9.5     user may amend as required for own model
#
load("c://source/mdvis.RData") 
poimd <- glm(numvisit ~ reform+ badh+age+educ+loginc,family=poisson, data=mdvis) 
mu <- fitted.values(poimd) 
p <- NULL 
avgp <- NULL
for (i in 0:20) { 
 p[[i+1]] <- exp(-mu)*(mu^i)/factorial(i) 
 avgp[i+1] <- mean(p[[i+1]])
} 
nCases <- dim(mdvis) 
n<- NULL 
propObs<- NULL 
probFit<- NULL 
yFitMean<- NULL 
for (i in 0:20) {                        #possible values for naffairs
bLos<- mdvis$numvisit==i                 #selector for naffairs=i 
n[i+1]<- sum(bLos)                       #number for naffairs=i 
propObs[i+1]<- n[i+1]/nCases[1]          #observed proportion for naffairs=i 
} 
Diff <- propObs*100 - avgp*100
data.frame(LOS=0:20, ObsProp=propObs*100, avgp*100, Diff)

plot(0:20, avgp, type="b", xlim=c(0,20), 
      main = "Observed vs Predicted Visits", 
      xlab = "Number Visits", ylab = "Probability of Visits")
lines(0:20, propObs, type = "b", pch = 2)
legend("topright", legend = c("Predicted Visits","Observed Visits"),
        lty = c(1,1), pch = c(1,2)) 




