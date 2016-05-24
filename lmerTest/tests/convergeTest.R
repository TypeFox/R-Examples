require(lmerTest)

###ASSERT ERROR

m.carrots <- lmer(Preference ~ sens1*sens2*Homesize*Age +
                     (1 + sens1 + sens2|Consumer) + (1|product), data = carrots)

step(m.carrots)

## model is not identifiable
tools::assertError(step(m.carrots, reduce.random = FALSE))

## check anova table with SAS
m.carrots <- lmer(Preference ~ sens1*sens2*Homesize*Age +
                    (1 +  sens2|Consumer) + (1|product), data = carrots)

an.carrots <- anova(m.carrots)

TOL <- 1e-5
stopifnot(all.equal(an.carrots$DenDF, c(10.6893, 13.9207, 94.8742, 94.8634,
                                        10.6881, 1003.08, 95.0939, 1003.09,
                                        95.0908, 94.8634, 1003.07, 1003.07, 
                                        1003.1, 95.0908, 1003.07), tol = TOL))








