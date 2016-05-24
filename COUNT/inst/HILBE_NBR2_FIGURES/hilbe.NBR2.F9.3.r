# hilbe.NBR2.F9.3.r
# Negative binomial regression postestimation graphic 
#    Standardized deviance vs fitted value, mu.  
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Figure 9.3     user may amend as required for own model
#
load("c://source/affairs.RData")
affnb2r <-  glm.nb(naffairs~ avgmarr + hapavg + vryhap + smerel 
     + vryrel + yrsmarr4 + yrsmarr5 + yrsmarr6, data=affairs) 
summary(affnb2r)
confint(affnbr2)
exp(coef(affnb2r))
exp(confint(affnb2r))
deviance <- residuals(affnb2, type="deviance")
dev <- sum(deviance*deviance)
pred <- predict(affnb2, se.fit=TRUE, type="response")
mu <- pred$fit
stdp <- pred$se.fit                 # Std error of prediction
variance <- mu 
h <- stdp * stdp*variance           # hat matrix diagonal
sdeviance <- rstandard(affnb2)      # Std deviance
plot(mu, sdeviance)

