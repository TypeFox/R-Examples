# fit some models
#
gpa.lm <- lm(gpa~satm+satv+act,gpa)
gpa.lma <- lm(gpa~ -1 + satm+satv+act,gpa)
#
# model comparison tests for 5 p-values in summary(gpa.lm)
#
anova(gpa.lma,gpa.lm)
