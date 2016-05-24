pet.lm <- lm(Rate ~ Group, petstress)
summary(Rate~Group,petstress,fun=favstats)
anova(pet.lm)
###hop:3-9
summary(pet.lm)
