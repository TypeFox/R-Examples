# hilbe.NBR2.F9.6ab.r
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# mu vs standardized deviance graphs
# Figure 9.6a and 9.6b     user may amend as required for own model
# Table 9.46  
#
# mu and Std deviance from ex4poie in T 9.44 - Poisson model
deviancep <- residuals(ex4poie, type="deviance")
devp <- sum(deviancep*deviancep)
predp <- predict(ex4poie, se.fit=TRUE, type="response")
mup <- predp$fit
mup <- 2*sqrt(mup)
stdpp <- predp$se.fit                 # Std error of prediction
variancep <- mup 
hp <- stdpp * stdpp*variancep           # hat matrix diagonal
sdeviancep <- rstandard(ex4poie)      # Std deviance
plot(mup, sdeviancep)
# mu and Std deviance from ex4nbe in T 9.45 - NB2 model
deviancenb <- residuals(ex4nbe, type="deviance")
devnb <- sum(deviancenb*deviancenb)
prednb <- predict(ex4nbe, se.fit=TRUE, type="response")
munb <- prednb$fit
munb <- 2*sqrt(munb)
stdpnb <- prednb$se.fit                 # Std error of prediction
variancenb <- munb 
hnb <- stdpnb * stdpnb*variancenb     
sdeviancenb <- rstandard(ex4nbe)       # Std dev 
plot(munb, sdeviancenb)                # Figure 9.6B

