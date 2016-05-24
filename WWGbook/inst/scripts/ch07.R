##############################################
#                                            #
#                 Chapter 7                  #
#                                            #
##############################################

 
library(WWGbook)

head(veneer)
attach(veneer)

age.f <- factor(age)
time.f <- factor(time)
tooth.f <- factor(tooth)

veneer <- data.frame(veneer,age.f, time.f, tooth.f)


###############################
#      Descriptive Analyses   #
###############################

sort(age.f)

#### Figure 7.2: Raw GCF values for each tooth vs. time, by patient. Panels are ordered by patient age.

library(lattice)  # Load the library for trellis graphics.
#trellis.device(color=F) # Make sure color is turned off.

# Load the nlme library, which is required for the 
# plots below as well as for subsequent models.
library(nlme)

veneer.g1 <- groupedData(gcf ~ time | tooth.f, 
outer = ~ age.f, data = veneer)
plot(veneer.g1, display = "tooth", outer = TRUE, aspect = 2, key = FALSE, xlab = "time", ylab = "GCF", main="age" ,
layout=c(4,3)) 

###############################
#      Model Analyses         #
###############################


library(nlme)

model7.1.fit <- lme(gcf ~ time + base.gcf + cda + age + time:base.gcf + time:cda + time:age, 
random = list(patient = ~ time, tooth = ~1), 
data = veneer, method = "REML")

summary(model7.1.fit)

intervals(model7.1.fit)

random.effects(model7.1.fit)

model7.1A.fit <- lme(gcf ~ time + base.gcf + cda + age + time:base.gcf + time:cda + time:age, 
random = list(patient = ~ time), data = veneer, method = "REML")

anova(model7.1.fit, model7.1A.fit)

model7.2A.fit <- lme(gcf ~ time + base.gcf + cda + age + time:base.gcf + time:cda + time:age, random = list(patient = ~ time, tooth = ~1), corr=corCompSymm(0.5, form=~1 | patient/tooth), weights = varIdent(form = ~1 | time), data = veneer, method = "REML")

summary(model7.2A.fit)

intervals(model7.2A.fit)

model7.2B.fit <- lme(gcf ~ time + base.gcf + cda + age + time:base.gcf + time:cda + time:age, random = list(patient = ~ time, tooth = ~1), corr=corCompSymm(0.5, form=~1 | patient/tooth), data = veneer, method = "REML")

summary(model7.2B.fit)

intervals(model7.2B.fit)
## Error in intervals.lme(model7.2B.fit)

model7.2C.fit <- lme(gcf ~ time + base.gcf + cda + age + time:base.gcf + time:cda + time:age, random = list(patient = ~ time, tooth = ~1), weights = varIdent(form = ~1 | time), data = veneer, method = "REML")

summary(model7.2C.fit)
intervals(model7.2C.fit)

anova(model7.2C.fit, model7.1.fit)

model7.1.ml.fit <- lme(gcf ~ time + base.gcf + cda + age + time:base.gcf + time:cda + time:age, 
random = list(patient = ~ time, tooth = ~1), data = veneer, method = "ML")

model7.3.ml.fit <- lme(gcf ~ time + base.gcf + cda + age, random = list(patient = ~ time, tooth = ~1), data = veneer, method = "ML")

anova(model7.1.ml.fit, model7.3.ml.fit)

model7.3.fit <- lme(gcf ~ time + base.gcf + cda + age, random = list(patient = ~ time, tooth = ~1), data = veneer, method = "REML")

summary(model7.3.fit)

############################
#  Diagnostics Analyses    #
############################

###Figure 7.4  
library(lattice)
trellis.device(color=F)
plot(model7.3.fit, resid(., type="p") ~ fitted(.),xlab="Predicted value", ylab="Residual",   abline=0)

###Figure 7.5
qqnorm(model7.3.fit, ~resid(.)) #layout=c(2,1), aspect = 2, id = 0.05) 

###Figure 7.6

qqnorm(model7.3.fit, ~ranef(.,level=1),layout=c(2,1), aspect = 2, id = 0.05) 


detach(veneer)

