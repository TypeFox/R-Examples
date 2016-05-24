## ----LIB, echo=FALSE, message=FALSE, warning=FALSE-----------------------
library(ANOM)
library(multcomp)
library(ggplot2)

## ----DATAHG, message=FALSE, warning=FALSE--------------------------------
library(ANOM)
data(hemoglobin)

## ----PSEUDOONEWAY, fig.cap="ANOM decision chart for the hemoglobin data based one a pseudo-one-way analysis.", fig.height=5.5----
hemoglobin$the <- as.factor(abbreviate(hemoglobin$therapy))
hemoglobin$td <- with(hemoglobin, the:drug)
hemodel <- lm(level ~ td, hemoglobin)
he <- glht(hemodel, mcp(td="GrandMean"), alternative="two.sided")
ANOM(he, xlabel="Treatment", ylabel="Hemoglobin Level")

## ----TWOWAYIA------------------------------------------------------------
hemodel2 <- lm(level ~ drug * therapy, hemoglobin)
anova(hemodel2)

## ----TWOWAY, fig.cap="ANOM decision chart for the hemoglobin data based on a two-way analysis.", fig.height=5.5----
hemodel3 <- lm(level ~ drug + therapy, hemoglobin)
he3 <- glht(hemodel3, mcp(drug="GrandMean"), alternative="two.sided")
ANOM(he3, xlabel="Drug", ylabel="Hemoglobin Level")

## ----INSECTSData---------------------------------------------------------
data(InsectSprays)
InsectSprays$block <- as.factor(rep(1:6, each=2))

## ----INSECTS-------------------------------------------------------------
insmodel1 <- glm(count ~ spray + block, data=InsectSprays,
                 family=quasipoisson(link="log"))
summary(insmodel1)$dispersion

## ----INSECTS2------------------------------------------------------------
insmodel2 <- glm(count ~ spray + block, data=InsectSprays,
                 family=poisson(link="log"))
anova(insmodel2, test="Chisq")

## ----INSECTS2a, fig.cap="ANOM decision chart for the insect spray data based on a Poisson GLM.", fig.height=5.5----
ins <- glht(insmodel2, mcp(spray="GrandMean"))
ANOM(ins)

## ----DATAES, message=FALSE, warning=FALSE--------------------------------
library(nlme)
data(ergoStool)

## ----NLME, fig.cap="ANOM decision chart for the ergonomic stool data based on a linear mixed-effects model.", fig.height=5.5----
library(nlme)
esmodel1 <- lme(effort ~ Type, random=~1|Subject, data=ergoStool)
es1 <- glht(esmodel1, mcp(Type="GrandMean"), alternative="two.sided")
ANOM(es1, xlabel="Stool Type", ylabel="Exertion (Borg Scale)")

## ----LME4, eval=FALSE----------------------------------------------------
#  library(lme4)
#  esmodel2 <- lmer(effort ~ Type + (1|Subject), data=ergoStool)
#  es2 <- glht(esmodel2, mcp(Type="GrandMean"), alternative="two.sided")
#  ANOM(es2, xlabel="Stool Type", ylabel="Exertion (Borg Scale)")

## ----MIXIGNORE, fig.cap="ANOM decision chart for the ergonomic stool data based on a standard linear model.", fig.height=5.5----
esmodel3 <- lm(effort ~ Type, ergoStool)
es3 <- glht(esmodel3, mcp(Type="GrandMean"), alternative="two.sided")
ANOM(es3, xlabel="Stool Type", ylabel="Exertion (Borg Scale)")

