##############################################
#                                            #
#                 Chapter 3                  #
#                                            #
##############################################

 
library(WWGbook)

head(ratpup)

###############################
#      Descriptive Analyses   #
###############################

#Figure 3.1:Boxplots to compare the distributions of birth weights 
####for each treatment by sex combination graphically

library(lattice)  # trellis graphics
library(grid)

bwplot(weight ~ sex|treatment, data=ratpup,aspect = 2, ylab="Birth Weights", 
    xlab="SEX",main = "Boxplots of birth weights for levels of treatment by sex")

#Figure 3.2:Boxplots to illustrate the relationship between birth weight and litter size. 
####Each panel shows the distributions of birth weights
####for all litters ranked by size, within a given level of treatment and sex. 

attach(ratpup)
ranklit <- litsize+0.01*litter
sort(ranklit)
ranklit
ranklit <- factor(ranklit)

levels(ranklit) <- as.character(1:27)
ranklit
bwplot(weight ~ranklit | treatment*sex, data=ratpup,aspect = 2, ylab="Birth Weights", 
    xlab="" )

###recode the variable SEX into SEX1
ratpup$sex1[sex == "Female"] <- 1
ratpup$sex1[sex == "Male"] <- 0

 
#########################################
#      Models Fitted by using lme       #
#########################################

library(nlme)

# Model 3.1.
model1.fit <- lme(weight ~ treatment + sex1 + litsize + treatment*sex1, 
    random = ~1 | litter, ratpup, method = "REML")
summary(model1.fit)
anova(model1.fit)

# Display the random effects (EBLUPs) from the model.
random.effects(model1.fit)

# Model 3.1A.
model1a.fit <- gls(weight ~ treatment + sex1 + litsize + treatment*sex1, data = ratpup)
anova(model1.fit, model1a.fit)  #Test Hypothesis 1.

# Model 3.2A.
model2a.fit <- lme(weight ~ treatment + sex1 + litsize + treatment*sex1, random = ~1 | litter, ratpup, method = "REML", weights = varIdent(form = ~1 | treatment))

# Test Hypothesis 3.2.
anova(model1.fit, model2a.fit)

ratpup$trtgrp[treatment == "Control"] <- 1
ratpup$trtgrp[treatment == "Low" | treatment == "High"] <- 2

model2b.fit <- lme(weight ~ treatment + sex1 + litsize + treatment*sex1, random = ~1 | litter, ratpup, method = "REML", weights = varIdent(form = ~1 | trtgrp))

# Test Hypothesis 3.3.
anova(model2a.fit, model2b.fit) 

# Test Hypothesis 3.4.
anova(model1.fit, model2b.fit) 

summary(model2b.fit) 

# Test Hypothesis 3.5.
anova(model2b.fit)  

model3.ml.fit <- lme(weight ~ treatment + sex1 + litsize, random = ~1 | litter, ratpup, method = "ML", weights = varIdent(form = ~1 | trtgrp))
model3a.ml.fit <- lme(weight ~ sex1 + litsize, random = ~1 | litter, ratpup, method = "ML", weights = varIdent(form = ~1 | trtgrp))

# Test Hypothesis 3.6.
anova(model3.ml.fit, model3a.ml.fit)

# Model 3.3: Final Model.
model3.reml.fit <- lme(weight ~ sex1 + litsize + treatment, random = ~1 | litter, ratpup, method = "REML", 
weights = varIdent(form = ~1 | trtgrp))
summary(model3.reml.fit)
anova(model3.reml.fit)


intervals(model3.reml.fit)

getVarCov(model3.reml.fit, individual="1", type="marginal")


############################
#  Model Diagnostics       #
############################


library(lattice)
trellis.device(color=F)

###Figure 3.4
res <- resid(model3.reml.fit)
ratred <- data.frame(ratpup, res)
ratred
histogram( ~res | factor(trtgrp),data=ratred, layout=c(2,1), aspect = 2 , xlab = "Residual") 

###Figure 3.5
qqnorm(model3.reml.fit, ~resid(.) | factor(trtgrp), layout=c(2,1), aspect = 2, id = 0.05) 

###Figure  3.6 
plot(model3.reml.fit, resid(., type="p") ~ fitted(.) | factor(trtgrp), layout=c(2,1), aspect=2, abline=0)

###Figure 3.7
bwplot(resid(model3.reml.fit) ~ranklit, data=model3.reml.fit, ylab="residual", xlab="litter")

detach(ratpup)

##########################################
#      Models Fitted by using lmer       #
##########################################

library(WWGbook)

library(Matrix)
 library(lme4) 
attach(ratpup)
ratpup$sex2[sex == "Female"] <- 1 
ratpup$sex2[sex == "Male"] <- 0 
# Model 3.1. 
model3.1.fit <- lmer(weight ~ treatment + sex2 + litsize + treatment*sex2 + (1 | litter), ratpup, method = "REML") 
summary(model3.1.fit) 
anova(model3.1.fit) 
ranef(model3.1.fit)$litter 

detach(ratpup)

