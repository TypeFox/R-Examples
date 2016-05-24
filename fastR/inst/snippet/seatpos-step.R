# trace=0 turns off intermediate reporting
seatpos.lmstep<-step(seatpos.lm1, trace=0)  
###hop:3-9
summary(seatpos.lmstep)
vif(seatpos.lmstep)
anova(seatpos.lm1,seatpos.lmstep)
