# Wong3-1  Political views, support for women to work and national welfare spending

library(vcdExtra)
# Data from Wong, R. (2010), Association Models, Los Angeles: Sage, Number 07-164
#      Table 3.1, from the General Social Survey, 2006.
# Questions:
#   polviews:  Think of self as liberal or conservative
#   fefam:     Better for men to work and women to tend home
#   natfare:   National welfare spending: too little, about right, too much

# Table 3.1
Freq<-c( 9, 5, 5, 1,    1, 6, 5, 1,    2, 2, 2, 1,
        17,13, 7, 4,   13,22, 9, 1,    7,13, 6, 2,
         8,14, 6, 0,   10,29,10, 0,    5,14, 6, 2,
        20,38,24, 8,   23,72,34,10,   17,67,36,12,
         4,21,12, 4,    7,30, 9, 1,    9,19,14, 2,
         2, 9, 8, 3,    1,16,19, 2,   11,28,28,11,
         0, 1, 5, 0,    2, 3, 3, 2,    2, 7, 6, 6)

polviews<-gl(7,4*3)
fefam<-gl(4,1,length=7*4*3)
natfare<-gl(3,4,length=7*4*3)

long.vnames <- list(set_varnames = c(polviews="Political views", fefam="Females should tend home",
                                     natfare="National welfare spending"))
long.lnames <- list(polviews = c("Lib++", "2", "3", "Moderate", "5", "6", "Cons++"),
                    fefam = c("Dis+", "Dis", "Agr", "Agr+"),
                    natfare = c("--", "OK", "++")
                    )


############################################
Wong31 <- data.frame(Freq, polviews, fefam, natfare)

Wong31.xtab <- xtabs(Freq ~ polviews+fefam+natfare, data=Wong31)
dimnames(Wong31.xtab) <- long.lnames

# Quick look at all pairwise associations
pairs(Wong31.xtab, gp=shading_Friendly, diag_panel=pairs_diagonal_mosaic)

############################################
# Model 1 - Independence Model
Wong31.O<-gnm(Freq~polviews+fefam+natfare,family=poisson, data=Wong31)
summary(Wong31.O)

mosaic(Wong31.O, main=paste("Independence model",modFit(Wong31.O)),
	labeling_args=long.vnames, set_labels=long.lnames,
	split_vertical=c(TRUE, FALSE, FALSE),
	labeling=labeling_residuals, suppress=2, gp=shading_Friendly)


############################################
# NB: add1 doesn't work with gnm() objects.  Re-fit using glm()
Wong31.O<-glm(Freq~polviews+fefam+natfare,family=poisson, data=Wong31)

# consider all two-way terms
add1(Wong31.O, ~.+(polviews + fefam + natfare)^2, test="Chisq")

# same result with MASS::addterm
addterm(Wong31.O, ~.+(polviews + fefam + natfare)^2, test="Chisq")

# or, start with saturated model and drop terms
Wong31.sat<-glm(Freq~polviews*fefam*natfare, family=poisson, data=Wong31)
drop1(Wong31.sat, test="Chisq")


############################################
# Model 2 - Full Two-way Interaction

Wong31.twoway <- update(Wong31.O, ~ .^2)
summary(Wong31.twoway)

mosaic(Wong31.twoway, main=paste("All two-way model", modFit(Wong31.twoway)),
	labeling_args=long.vnames, set_labels=long.lnames,
	split_vertical=c(TRUE, FALSE, FALSE),
	labeling=labeling_residuals, suppress=1, gp=shading_Friendly)


############################################
# Model 3 - Conditional Independence on polviews
Wong31.cond1 <- glm(Freq~polviews * (fefam + natfare), family=poisson)
summary(Wong31.cond1)

mosaic(Wong31.cond1, main=paste("Cond1: ~P * (F+N)", modFit(Wong31.cond1)),
	labeling_args=long.vnames, set_labels=long.lnames,
	split_vertical=c(TRUE, FALSE, FALSE),
	labeling=labeling_residuals, suppress=1, gp=shading_Friendly)

############################################
# Model 4 - Conditional Independence on fefam
Wong31.cond2 <- glm(Freq~polviews*fefam + fefam*natfare,family=poisson)
summary(Wong31.cond2)

mosaic(Wong31.cond2, main=paste("Cond2: ~F * (P+N)", modFit(Wong31.cond2)),
	labeling_args=long.vnames, set_labels=long.lnames,
	split_vertical=c(TRUE, FALSE, FALSE),
	labeling=labeling_residuals, suppress=1, gp=shading_Friendly)

############################################
# Model 5 - Conditional Independence on natfare
Wong31.cond3<-glm(Freq~fefam*natfare+polviews*natfare,family=poisson)
summary(Wong31.cond3)

mosaic(Wong31.cond3, main=paste("Cond2: ~N * (F+N)", modFit(Wong31.cond3)),
	labeling_args=long.vnames, set_labels=long.lnames,
	split_vertical=c(TRUE, FALSE, FALSE),
	labeling=labeling_residuals, suppress=1, gp=shading_Friendly)

anova(Wong31.O,  Wong31.cond3, Wong31.cond2, Wong31.cond1,Wong31.twoway, Wong31.sat)
