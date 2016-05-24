dome.lm1 <- lm(Dist~Velocity+Angle+BallWt+BallDia+Cond, data=domedata)
step(dome.lm1,direction="both")
dome.lm2 <- lm(Dist~Velocity+Cond,data=domedata)
###hop:3-9
summary(dome.lm2)
anova(dome.lm2,dome.lm1)
dome.lm3 <- lm(Dist~Velocity+Cond+BallDia,data=domedata)
anova(dome.lm2,dome.lm3)
