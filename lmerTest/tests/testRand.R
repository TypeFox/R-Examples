##################################################################
## check elimination of random effects in rand and step functions
##################################################################

require(lmerTest)

## check for random coefficient models
modelCarrots <- lmer(Preference ~  sens2*sens1*Homesize*Age
                                    + (1 | product) + (1 + sens1 + sens2 | Consumer),
                                     data=carrots)

## the results for the rand function differ from the step
## because of the update function - changing the contrasts to contr.SAS
rnd <- rand(modelCarrots)



stopifnot(all.equal(rnd$rand.table[,"p.value"], c(0.0000289113, 0.4852719557,
                                                  0.0186101913)))

stp <- step(modelCarrots)

stp




modelCarrots_red1 <- lmer(Preference ~ sens2*sens1*Homesize*Age
                                     + (1 | product) + (1   + sens2 | Consumer),
                                     data=carrots)

modelCarrots_red2 <- lmer(Preference ~
                            sens2*sens1*Homesize*Age
                          + (1 | product) + (1 | Consumer),
                          data=carrots)


modelCarrots_redp <- lmer(Preference ~
                            sens2*sens1*Homesize*Age
                           + (1   + sens2| Consumer),
                          data=carrots)

res.lrt.s1 <- anova(modelCarrots, modelCarrots_red1, refit = FALSE)
res.lrt.s2 <- anova(modelCarrots_red1, modelCarrots_red2, refit = FALSE)
res.lrt.sp <- anova(modelCarrots_red1, modelCarrots_redp, refit = FALSE)


stopifnot(all.equal(res.lrt.s1[2,"Pr(>Chisq)"], stp$rand.table[1,"p.value"]),
          all.equal(res.lrt.s2[2,"Pr(>Chisq)"], stp$rand.table[3,"p.value"]),
          all.equal(res.lrt.sp[2,"Pr(>Chisq)"], stp$rand.table[2,"p.value"]))


