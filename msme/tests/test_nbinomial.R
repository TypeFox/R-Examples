# nbinomial.r    msme package
#                Joseph M Hilbe and Andrew P Robinson
#                Methods of Statistical Model Estimation
#                Chapman & Hall/CRC  2013
# Examples of use for the nbinomial.r function for maximum likelihood 
# NB2 and heterogeneous NB2 regression. The default "family" is 
# "nb2", which estimates parameters and the dispersion parameter using 
# a direct relationship between the dispersion and Poisson variance. 
# The default NB variance is mu + alpha*mu^2. Use the "negBinomial" family
# to produce the inverted dispersion, which is given by R's glm.nb function

library(msme)

# library(msme, lib.loc="lib")

data(medpar)

# TRADITIONAL NB REGRESSION WITH ALPHA
mynb1 <- nbinomial(los ~ hmo + white, data=medpar)
summary(mynb1)

# DISPERSION STATISTIC A
mynb1$dispersion

# INCIDENCE RATE RATIOS
exp(mynb1$coefficients)

# IRR SEs USING DELTA METHOD
exp(mynb1$coefficients)*mynb1$se.beta.hat

#IRR CONFIDENCE INTERVALS
mynb1$coefficients - 1.96*mynb1$se.beta.hat
mynb1$coefficients + 1.96*mynb1$se.beta.hat

# TRADITIONAL NB -- SUMMARY INCLUDED IN FUNCTION CALL
summary(mynb1_5 <- nbinomial(los ~ hmo + white, data=medpar))

# TRADITIONAL NB -- SHOWING ALL OPTIONS
mynb2 <- nbinomial(los ~ hmo + white,
                    formula2 = ~ 1,
                    data = medpar,
                    family = "nb2",
                    mean.link = "log",
                    scale.link = "inverse_s")
summary(mynb2)

# R GLM.NB - LIKE INVERTED DISPERSION BASED M
mynb3 <- nbinomial(los ~ hmo + white,
                    formula2 = ~ 1,
                    data = medpar,
                    family = "negBinomial",
                    mean.link = "log",
                    scale.link = "inverse_s")
summary(mynb3)

# R GLM.NB-TYPE INVERTED DISPERSON --THETA ; WITH DEFAULTS
mynb4 <- nbinomial(los ~ hmo + white, family="negBinomial", data =medpar)
summary(mynb4)

# HETEROGENEOUS NB; DISPERSION PARAMETERIZED
mynb5 <- nbinomial(los ~ hmo + white,
                    formula2 = ~ hmo + white,
                    data = medpar,
                    family = "negBinomial",
                    mean.link = "log",
                    scale.link = "log_s")
summary(mynb5)

# SAVED STATISTICS FROM NBINOMIAL2.R


# USE OF "nbh" for glm.nb-like dipersion
# nbinomial now has "negBinomial" as the family for glm.nb


# TRADITIONAL NB2 REGRESSION. USING DEFAULT OPTIONS
tnb1 <- nbinomial(los ~ hmo + white, data=medpar)
summary(tnb1)


# TRADITIONAL NB2 REGRESSION. PROVIDING EXPLICIT OPTIONS
tnb2 <- nbinomial(los ~ hmo + white,
                    formula2 = ~ 1,
                    data = medpar,
                    family = "nb2",
                    mean.link = "log",
                    scale.link = "inverse_s")

summary(tnb2)

# PEARSON CHI2 STATISTIC
tnb1$pearson



