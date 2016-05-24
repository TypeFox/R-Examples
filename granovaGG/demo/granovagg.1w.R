data(arousal)
#Drug A
granovagg.1w(arousal[,1:2], h.rng = 1.6, v.rng = 0.5)

###

data(anorexia)
wt.gain <- anorexia[, 3] - anorexia[, 2]
granovagg.1w(wt.gain, group = anorexia[, 1])

###

data(poison)
##Note violation of constant variance across groups in following graphic.
granovagg.1w(poison$SurvTime, group = poison$Group, ylab = "Survival Time")
##RateSurvTime = SurvTime^-1
granovagg.1w(poison$RateSurvTime, group = poison$Group,
ylab = "Survival Rate = Inverse of Survival Time")

##Nonparametric version: RateSurvTime ranked and rescaled
##to be comparable to RateSurvTime;
##note labels as well as residual (rug) plot below.
granovagg.1w(poison$RankRateSurvTime, group = poison$Group,
ylab = "Ranked and Centered Survival Rates",
main = "One-way ANOVA display, poison data (ignoring 2-way set-up)",
res = TRUE)

###

data(chickwts)
?chickwts # An explanation of the chickwts dataset
with(chickwts, granovagg.1w(weight, group = feed)) # Modeling weight as explained by feed type
