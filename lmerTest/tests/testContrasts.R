## TODO: load data tree here

require(lmerTest)

load(system.file("testdata", "tree.RData", package="lmerTest"))

modelCarrots.treat <- lme4::lmer(Preference ~
                                   sens2*sens1*Homesize*Age
                                 + (1 | product) + 
                                   (1   + sens1 + sens2 | Consumer),
                                 data=carrots, 
                                 contrasts = list(Homesize = "contr.treatment", 
                                                  Age = "contr.treatment"))

modelCarrots.sas <- lme4::lmer(Preference ~
                                   sens2*sens1*Homesize*Age
                                 + (1 | product) + 
                                 (1  + sens1 + sens2 | Consumer),
                                 data=carrots, 
                              contrasts = list(Homesize = "contr.SAS", 
                                               Age = "contr.SAS"))

## here an error produces
## NO MORE in lme4 1.1-8
# tools::assertError(stopifnot(all.equal(logLik(modelCarrots.treat), 
#                                        logLik(modelCarrots.sas))))
# 
# tools::assertError(stopifnot(all.equal(VarCorr(modelCarrots.treat), 
#                                        VarCorr(modelCarrots.sas), tol = 1e-5)))


                              
modelHam.sas <- lmer(Informed.liking ~
                 Product*Information*Gender*Age
                 + (1 | Consumer)
                 + (1 | Product:Consumer)
                 + (1 | Information:Consumer), data=ham, 
                 contrasts = list(Product = "contr.SAS", 
                                  Information = "contr.SAS", 
                                  Gender = "contr.SAS"))

modelHam.treat <- lmer(Informed.liking ~
                       Product*Information*Gender*Age
                     + (1 | Consumer)
                     + (1 | Product:Consumer)
                     + (1 | Information:Consumer), data=ham, 
                     contrasts = list(Product = "contr.treatment", 
                                      Information = "contr.treatment", 
                                      Gender = "contr.treatment"))

stopifnot(all.equal(logLik(modelHam.sas), logLik(modelHam.treat)))
stopifnot(all.equal(VarCorr(modelHam.sas), VarCorr(modelHam.treat), 
                    tol = 1e-3))

## check that lsmeans is the same whether the contrasts for the models are differenr
lmer4 <- lmer(increase ~ treat + (1|block), data = tree,  
              contrasts = list(treat = "contr.treatment"))

lmer5 <- lmer(increase ~ treat+ (1|block), data = tree,  
              contrasts = list(treat = "contr.SAS"))

all.equal(lsmeans(lmer4), lsmeans(lmer5), tol = 1e-3)

