# Wong2-3  Political views and support for women to work

library(vcdExtra)
# Data from Wong, R. (2010), Association Models, Los Angeles: Sage, Number 07-164
#      Table 2.3A, from the General Social Survey, 1998-2000.
# Questions:
#   polviews:  Think of self as liberal or conservative
#   fefam:     Better for men to work and women to tend home

Freq<-c(39, 50, 18,  4,
       140,178, 85, 23,
       108,195, 97, 23,
       238,598,363,111,
        78,250,150, 55,
        50,200,208, 74,
         8, 29, 46, 21)
polviews<- gl(7,4)
fefam <- gl(4,1,length=28)

# create better labels for levels in mosaic()
## but, this screws up parameter names: use set_labels instead
#polviews <- factor(polviews, labels=c("Lib++", "2", "3", "Moderate", "5", "6", "Cons++"))
#fefam <- factor(fefam, labels=c("Dis+", "Dis", "Agr", "Agr+"))
#
long.vnames <- list(set_varnames = c(polviews="Political views", fefam="Females should tend home"))
long.lnames <- list(polviews = c("Lib++", "2", "3", "Moderate", "5", "6", "Cons++"),
                    fefam = c("Dis+", "Dis", "Agr", "Agr+"))


# numeric versions for U, R, C, RC models
Rscore<-as.numeric(polviews)
Cscore<-as.numeric(fefam)

# make a data frame
Wong23 <- data.frame(Freq, polviews, fefam, Rscore, Cscore)

####################################
# do a correspondence analysis first, to see the row/category relations
Wong23.xtab <- xtabs(Freq ~ polviews+fefam, data=Wong23)
dimnames(Wong23.xtab) <- long.lnames
library(ca)

plot(ca(Wong23.xtab))
title(main="Political views and support for women to work", 
      xlab="Dim 1 (90.8%)", ylab="Dim 2 (8.5%)")


####################################


## OK now, gives warning
Wong23.O <- gnm(Freq~polviews+fefam, family=poisson, data=Wong23)
 
# OK, with formula
mosaic(Wong23.O, main=paste("Independence model", modFit(Wong23.O)), 
	formula=~polviews+fefam,
	labeling_args=long.vnames, set_labels=long.lnames)

####################################
# Uniform association model

Wong23.U<-gnm(Freq~polviews+fefam+Rscore:Cscore,family=poisson,tolerance = 1e-12, data=Wong23)
anova(Wong23.U)
# OK, w/o formula
mosaic(Wong23.U, main=paste("Uniform association", modFit(Wong23.U)),
	formula=~polviews+fefam,
	labeling_args=long.vnames, set_labels=long.lnames)

# display standardized residuals
mosaic(Wong23.U, formula=~polviews+fefam, main=paste("Uniform association", modFit(Wong23.U)),
	labeling_args=long.vnames, set_labels=long.lnames, residuals_type="rstandard",
	labeling=labeling_residuals, suppress=1)

# compare with gp=shading_Friendly)
mosaic(Wong23.U, formula=~polviews+fefam, main=paste("Uniform association", modFit(Wong23.U)),
	labeling_args=long.vnames, set_labels=long.lnames, residuals_type="rstandard",
	labeling=labeling_residuals, suppress=1, gp=shading_Friendly)

####################################

# Model B - R: Row Effects
Wong23.R<-gnm(Freq~polviews+fefam+Cscore:polviews,family=poisson, data=Wong23)
anova(Wong23.R)
mosaic(Wong23.R, formula=~polviews+fefam, main="Row effects model",
	labeling_args=long.vnames, set_labels=long.lnames, residuals_type="rstandard",
	labeling=labeling_residuals, suppress=1, gp=shading_Friendly)

#####################################
# Model C - C: Column Effect
Wong23.C<-gnm(Freq~polviews+fefam+Rscore:fefam,family=poisson, data=Wong23)
anova(Wong23.C)
mosaic(Wong23.C, formula=~polviews+fefam, main="Column effects model",
	labeling_args=long.vnames, set_labels=long.lnames, residuals_type="rstandard",
	labeling=labeling_residuals, suppress=1, gp=shading_Friendly)
	
#####################################
# Model D - R+C: Row and Column Effect
oldopt <- options(contrasts = c(factor="contr.treatment", ordered="contr.treatment"))
Wong23.RplusC<-gnm(Freq~polviews+fefam+Rscore:Cscore+Cscore:polviews+Rscore:fefam,
         constrain=c(17,20),constrainTo=c(0,0),
         family=poisson,tolerance = 1e-12)
anova(Wong23.RplusC)
mosaic(Wong23.RplusC, formula=~polviews+fefam, main="Column effects model",
	labeling_args=long.vnames, set_labels=long.lnames, residuals_type="rstandard",
	labeling=labeling_residuals, suppress=1, gp=shading_Friendly)

options(oldopt)

#####################################
# Model E - RC: RC(1) model
Wong23.RC1<-gnm(Freq~polviews+fefam+Mult(1,polviews,fefam),
         family=poisson,tolerance = 1e-12)

mosaic(Wong23.RC1, formula=~polviews+fefam, main="RC(1) model",
	labeling_args=long.vnames, set_labels=long.lnames, residuals_type="rstandard",
	labeling=labeling_residuals, suppress=1, gp=shading_Friendly)

#####################################
# LRstats the collection of models
models <- list(Indep=Wong23.O, Uniform=Wong23.U, RowEff=Wong23.R, ColEff=Wong23.C, RplusC=Wong23.RplusC, RC1=Wong23.RC1)
res <- lapply(models, residuals)
boxplot(as.data.frame(res), main="Residuals from various models")

aic <- t(as.data.frame(lapply(models, extractAIC)))
colnames(aic) <- c("df", "AIC")
aic
# sort by df
aic <- aic[order(aic[,1]),]
plot(aic, type = "b", main="AIC plot")
text(aic, labels=rownames(aic), pos=c(4,1,3,1,3,1))

######################################
# compare models;  they are not nested, so only some Chisq tests make sense
anova(Wong23.O, Wong23.U, Wong23.R, Wong23.C, Wong23.RplusC, Wong23.RC1)

anova(Wong23.O, Wong23.U, Wong23.R, Wong23.RplusC, test="Chisq")
anova(Wong23.O, Wong23.U, Wong23.C, Wong23.RplusC, test="Chisq")

