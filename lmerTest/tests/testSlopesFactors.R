require(lmerTest)

load(system.file("testdata","bugSummaryData.RData", package="lmerTest"))

lmer3 <- lmer(cog_Mid ~ (1|item) + (1 + vowel|speaker) + 
                Language*vowel*sex, 
              data=data1.frame, REML = FALSE, na.action = na.omit)

## no more errors if the C code from lme4 is used
s <- summary(lmer3)

lmer4 <- lmer(cog_Mid ~ (1|item) + (vowel - 1|speaker) + Language*vowel*sex, 
              data = data1.frame)

## anova does not work in this case A is not positiv definite  - FIXED
## was no a positive definit because of Cv_to_Sv function...?
## SAS seems to work
## still some disagreemants in estimation of deviance, residuals
an <- anova(lmer4)


TOL <- 1e-3
stopifnot(all.equal(an$DenDF, c(61.0403, 36.5587, 67.0312, 68.505,
                                134.152, 81.2094, 215.628), tol = TOL))

stopifnot(all.equal(an$F.value, c(40.6522, 6.2784, 46.0567, 
                                  4.2176, 2.7903, 2.5043, 0.2995), tol = TOL))
