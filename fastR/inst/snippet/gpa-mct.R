# fit some models
#
gpa.lm <- lm(gpa~satm+satv+act,gpa)
gpa.lma <- lm(gpa~ -1 + satm+satv+act,gpa)
gpa.lmb <- lm(gpa~satv+act,gpa)
gpa.lmc <- lm(gpa~satm+act,gpa)
gpa.lmd <- lm(gpa~satm+satv,gpa)
gpa.lme <- lm(gpa~1,gpa)
#
# model comparison tests for 5 p-values in summary(gpa.lm)
#
anova(gpa.lma,gpa.lm)
anova(gpa.lmb,gpa.lm)
anova(gpa.lmc,gpa.lm)
anova(gpa.lmd,gpa.lm)
anova(gpa.lme,gpa.lm)
summary(gpa.lm)

#
# combined SAT verses subscore
#
gpa.lmf <- lm(gpa~I(satv+satm) + act,gpa)
anova(gpa.lmf,gpa.lm)
