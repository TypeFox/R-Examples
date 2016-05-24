# hilbe.NBR2.F14.1.r
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# FIGURE 14.1 
library(gamlss.mx)
load("c://source/medpar.RData")
rinb <- gamlssNP(los~ hmo +white+ type2 +type3, random=~1|provnum, 
      data=medpar, family=NBI, mixture="gq", K=20) 
summary(rinb)

summary(rinb$sigma.fv) 
m<-rinb$mu.fv       # fitted values for extended model
s<-rinb$sigma.fv    # sigma
presid <- (medpar$los-m)/sqrt(m+s*m*m) 
summary(presid)
hist(presid)





