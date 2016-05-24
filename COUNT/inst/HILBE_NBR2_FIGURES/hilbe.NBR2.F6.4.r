# hilbe.NBR2.F6.4.r
# Conditional effects Poisson plot with user specified mean values
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Table 6.17; Figure 6.4  
# User to amend data, model variables, and effects for graphing
load("c://source/rwm5yr.RData")
eppoi <- glm(docvis ~ outwork+age+female+married+edlevel2+edlevel3+ edlevel4,  family=poisson, data=rwm5yr)
rest  <- eppoi$coef[4]*mean(rwm5yr$female) + eppoi$coef[5]*mean(rwm5yr$married) +  
         eppoi$coef[6]*mean(rwm5yr$edlevel2) + eppoi$coef[7]*mean(rwm5yr$edlevel3) +  
         eppoi$coef[8]*mean(rwm5yr$edlevel4)
out0  <- eppoi$coef[1] + eppoi$coef[3]*rwm5yr$age + rest
out1  <- eppoi$coef[1] + eppoi$coef[2]*1 + eppoi$coef[3]*rwm5yr$age + rest
eout1 <- exp(out1)
eout0 <- exp(out0)
matplot(cbind(rwm5yr$age, rwm5yr$age), cbind(eout0, eout1), 
  pch=1:2, col=1:2, xlab='Count', ylab='Frequency?') 
matplot(cbind(rwm5yr$age, rwm5yr$age), cbind(eout0, eout1), type='l', 
  lty=1:2, col=1:2, xlab='Doctor visits', ylab='Frequency')
