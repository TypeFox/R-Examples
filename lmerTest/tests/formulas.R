require(lmerTest)

sens13 <- (carrots$sens1)^3
fm2 <- lmer(Preference ~ sens1 + I(sens1^2) + sens13 + Homesize +
              +             (1+sens2|Consumer), data=carrots)

an <- anova(fm2)

TOL <- 1e-3
stopifnot(all.equal(an[,"Pr(>F)"] , c(0.0388, 0.0002, 0.0814, 0.0202) , tol=TOL), 
          all.equal(round(an[,"DenDF"]) , c(1125, 1101, 1127, 101) , tol=TOL), 
          TRUE)

fm3 <- lmer(Preference ~ poly(sens1, 3)+ Homesize +
              +             (1+sens2|Consumer), data=carrots)

anova(fm3)


s.lmerTest <- summary(fm3)

s.lme4<- summary(fm3, ddf="lme4")

all.equal(coefficients(s.lmerTest)[, "t value"], 
          coefficients(s.lme4)[, "t value"])
