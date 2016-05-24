# UCBAdmissions data: Conditional independence via loglm and glm
library(vcd)
data("UCBAdmissions")

structable(Dept ~ Admit+Gender,UCBAdmissions)

## conditional independence in UCB admissions data
mod.1 <- loglm(~ Dept * (Gender + Admit), data=UCBAdmissions)
mod.1

# this is correct, except the Pearson residuals dont show that
# all the lack of fit is concentrated in Dept A
mosaic(mod.1, gp=shading_Friendly, labeling=labeling_residuals)


library(vcdExtra)

# using glm()
berkeley <- as.data.frame(UCBAdmissions)
mod.3 <- glm(Freq ~ Dept * (Gender+Admit), data=berkeley, family="poisson")
summary(mod.3)

# (BUG FIXED )the large residuals are all in Dept A
mosaic(mod.3, residuals_type="rstandard", labeling=labeling_residuals)

