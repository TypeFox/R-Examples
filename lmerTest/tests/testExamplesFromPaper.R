require(lmerTest)

modelTVbo <- lmer(Colourbalance ~ TVset*Picture
                  + (1 | Assessor) + (1 | Assessor:TVset)
                  + (1 | Assessor:Picture)
                  + (1 | Assessor:TVset:Picture)
                  + (1 | Repeat) + (1 | Repeat:TVset)
                  + (1 | Repeat:Picture) + (1|TVset:Picture:Repeat), data=TVbo)

stTV <- step(modelTVbo)

stTV

modelTVbo_red <- lmer(Colourbalance ~ TVset*Picture
                        + (1 | Assessor:TVset)
                        + (1 | Assessor:Picture)
                        +  (1 | Repeat:TVset)
                       , data=TVbo)

modelTVbo_redAT <- lmer(Colourbalance ~ TVset*Picture
                      + (1 | Assessor:Picture)
                      +  (1 | Repeat:TVset)
                      , data=TVbo)

modelTVbo_redAP <- lmer(Colourbalance ~ TVset*Picture
                      + (1 | Assessor:TVset)
                      +  (1 | Repeat:TVset)
                      , data=TVbo)

modelTVbo_redRT <- lmer(Colourbalance ~ TVset*Picture
                      + (1 | Assessor:TVset)
                      + (1 | Assessor:Picture)
                      , data=TVbo)

anova(modelTVbo_red, modelTVbo_redAT, refit = FALSE)

stopifnot(all.equal(anova(modelTVbo_red, 
                          modelTVbo_redAP, 
                          refit = FALSE)[2,"Pr(>Chisq)"],
                    stTV$rand.table["Assessor:Picture","p.value"]))

stopifnot(all.equal(anova(modelTVbo_red, 
                          modelTVbo_redAT, 
                          refit = FALSE)[2,"Pr(>Chisq)"],
                    stTV$rand.table["Assessor:TVset","p.value"]))

stopifnot(all.equal(anova(modelTVbo_red, 
                          modelTVbo_redRT, 
                          refit = FALSE)[2,"Pr(>Chisq)"],
                    stTV$rand.table["Repeat:TVset","p.value"]))


## for ham data
modelHam <- lmer(Informed.liking ~
                 Product*Information*Gender*Age
                 + (1 | Consumer)
                 + (1 | Product:Consumer)
                 + (1 | Information:Consumer), data=ham)

stHam <- step(modelHam)

stHam 

TOL <- 1e-5
stopifnot(all.equal(as.vector(stHam$lsmeans.table[1:4, "Estimate"]), 
                    as.vector(tapply(ham$Informed.liking, ham$Product, mean))
                    , tol = TOL,
                    check.attributes = FALSE, check.names = FALSE))

## the carrots example is tested in testRand.R function

