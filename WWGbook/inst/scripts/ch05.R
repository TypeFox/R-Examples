##############################################
#                                            #
#                 Chapter 5                  #
#                                            #
##############################################
library(WWGbook)
head(rat.brain)
attach(rat.brain)

###############################
#      Descriptive Analyses   #
###############################

region.f <- region
region.f[region == 1] <- 1
region.f[region == 2] <- 2
region.f[region == 3] <- 0
region.f <- factor(region.f)
treat <- treatment
treat[treatment == 1] <- 0
treat[treatment == 2] <- 1
treat <- factor(treat)
rat.brain <- data.frame(rat.brain, region.f, treat)

#######table report on page 221######

###overall ###
dim(rat.brain)  #obs
summary(rat.brain$activate) #mean
sqrt(var(rat.brain$activate)) #stdev

###by treatment group###
summary(treat)     #obs
by(rat.brain$activate,list(treat),summary) #mean
sqrt(by(rat.brain$activate,list(treat),var))   #stdev

###by region group###
summary(region.f) #obs
by(rat.brain$activate,list(region.f),summary) #mean
sqrt(by(rat.brain$activate,list(region.f),var))   #stdev

###by treat and region group###
table(treat,region.f)  #obs
by(rat.brain$activate,list(treat, region.f),summary)  #mean
sqrt(by(rat.brain$activate,list(treat, region.f),var)) #stdev

###Figure 5.1: Line graphs of activation for each animal by region within levels of treatment for the Rat Brain data.

library(lattice)  # Load the library for trellis graphics.
#trellis.device(color=F) # Make sure color is turned off.

# Load the nlme library, which is required for the 
# plots below as well as for subsequent models.
library(nlme)

rat.brain.g1 <- groupedData(activate ~ region | animal, 
outer = ~ treat, data = rat.brain)

plot(rat.brain.g1, display = "animal", outer = TRUE, aspect = 2, key = FALSE, xlab = "region", ylab = "mean activate", main="treatment" ) 

 
###############################
#      Models Fitted          #
###############################
 

# Model 5.1.
model5.1.fit <- lme(activate ~ region.f*treat, 
random = ~1 | animal, method = "REML", data = rat.brain)

summary(model5.1.fit)
anova(model5.1.fit)

# Model 5.2.
model5.2.fit <- update(model5.1.fit, random = ~treat | animal)

summary(model5.2.fit)
anova(model5.2.fit)

0.5*(1 - pchisq(26.1,1)) + 0.5*(1 - pchisq(26.1,2))


# Model 5.3.
model5.3.fit <- lme(activate ~ region.f*treat, 
random = ~treat | animal, 
weights = varIdent(form = ~1 | treat), 
data = rat.brain)

summary(model5.3.fit)

# Likelihood ratio test for Hypothesis 5.2.
anova(model5.2.fit, model5.3.fit)

getVarCov(model5.2.fit, individual = "R100797", type = "marginal")

############################
#  Model Diagnostics       #
############################

###Figure 5.3
qqnorm(model5.2.fit, ~resid(.)) #layout=c(2,1), aspect = 2, id = 0.05) 

###Figure 5.4  
fit <- predict(model5.2.fit)
res <- resid(model5.2.fit)
plotres.fit <-data.frame(rat.brain, fit, res)
xyplot(res~fit,data=plotres.fit,groups=treatment, xlab="Predicted value", ylab="Residual",abline=0)


detach(rat.brain)

