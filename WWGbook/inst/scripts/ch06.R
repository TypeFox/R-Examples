##############################################
#                                            #
#                 Chapter 6                  #
#                                            #
##############################################
library(WWGbook)
head(autism)
attach(autism)

sicdegp.f <- factor(sicdegp)
age.f <- factor(age)

# Add the new variables to a new data frame object.
autism.updated <- data.frame(autism, sicdegp.f, age.f)

###############################
#      Descriptive Analyses   #
###############################

# Number of Observations at each level of AGE
summary(age.f)
# Number of Observations at each level of AGE within each group
# defined by the SICDEGP factor
table(sicdegp.f, age.f)
# Overall summary for VSAE
summary(vsae)
# VSAE means at each AGE
tapply(vsae, age.f, mean, na.rm=TRUE)
# VSAE minimum values at each AGE
tapply(vsae, age.f, min, na.rm=TRUE)
# VSAE maximum values at each AGE
tapply(vsae, age.f, max, na.rm=TRUE)


#### Figure 6.1: Observed VSAE values plotted against age for children in each SICD group. 

library(lattice)  # Load the library for trellis graphics.
#trellis.device(color=F) # Make sure color is turned off.

# Load the nlme library, which is required for the 
# plots below as well as for subsequent models.
library(nlme)

autism.g1 <- groupedData(vsae ~ age | childid, 
outer = ~ sicdegp.f, data = autism.updated)

plot(autism.g1, display = "childid", outer = TRUE, aspect = 2, key = FALSE, xlab = "Age (Years)", ylab = "VSAE", 
main = "Individual Data by SICD Group") 

#### Figure 6.2: Mean Profiles of VSAE values for children in each SICD group. 

autism.g2 <- groupedData(vsae ~ age | sicdegp.f, 
order.groups = F, data = autism.updated)

plot(autism.g2, display = "sicdegp", aspect = 2, key = FALSE, 
xlab = "Age (Years)", ylab = "VSAE", 
main = "Mean Profiles by SICD Group")

###############################
#      Models Fitted          #
###############################

# Compute age.2 (AGE minus 2) and age.2sq (age.2 squared).
age.2 <- age - 2
age.2sq <- age.2*age.2

# Recode the SICDEGP factor for model fitting.
sicdegp2 <- sicdegp
sicdegp2[sicdegp == 3] <- 0
sicdegp2[sicdegp == 2] <- 2 
sicdegp2[sicdegp == 1] <- 1
sicdegp2.f <- factor(sicdegp2)

autism.updated <- subset(data.frame(autism, sicdegp2.f, age.2), !is.na(vsae))

library(nlme)

autism.grouped <- groupedData(vsae ~ age.2 | childid, data=autism.updated, order.groups = F)

model1.fit <- lme(vsae ~ age.2 + I(age.2^2) + sicdegp2.f +
age.2:sicdegp2.f + I(age.2^2):sicdegp2.f, 
random = ~ age.2 + I(age.2^2), method="REML", 
data = autism.grouped)
##Error: not converge

summary(model1.fit)

model2.fit <- lme(vsae ~ age.2 + I(age.2^2) + sicdegp2.f +
age.2:sicdegp2.f + I(age.2^2):sicdegp2.f, 
random = ~ age.2 + I(age.2^2) - 1, method="REML", 
data = autism.grouped)

summary(model2.fit)


model2a.fit <- update(model2.fit, random = ~ age.2 - 1)

h1.pvalue <- 0.5*(1-pchisq(83.9,1)) + 0.5*(1-pchisq(83.9,2))
h1.pvalue


model2.ml.fit <- update(model2.fit, method = "ML")

model3.ml.fit <- update(model2.ml.fit, 
fixed = ~ age.2 + I(age.2^2) + sicdegp2.f + age.2:sicdegp2.f)

anova(model2.ml.fit, model3.ml.fit)

model3.fit <- update(model2.fit, 
fixed = ~ age.2 + I(age.2^2) + sicdegp2.f + age.2:sicdegp2.f)

summary(model3.fit)

intervals(model3.fit)

getVarCov(model3.fit, individual="1", type="marginal")

############################
#  Model Diagnostics       #
############################

###Figure 6.4
curve(0.11*x^2 + 6.15*x + 13.46, 0, 11, xlab = "AGE minus 2", ylab = "Marginal Predicted VSAE", lty = 3, ylim=c(0,100), lwd = 2)

curve(0.11*x^2 + 2.65*x + 9.84, 0, 11, add=T, lty = 2, lwd = 2)

curve(0.11*x^2 + 2.08*x + 8.47, 0, 11, add=T, lty = 1, lwd = 2)

legend(locator(1),c("SICDEGP = 1", "SICDEGP = 2", 
"SICDEGP = 3"), lty = c(1,2,3), lwd = c(2,2,2))

###click the mouse button on the Figure 6.4 to put the legend on the place where you want


library(lattice)
trellis.device(color=F)

###Figure 6.5
plot(augPred(model3.fit, level = 0:1), layout=c(4,3,1),
xlab="AGE minus 2", ylab="Predicted VSAE",
key = list(lines = list(lty = c(1,2), col = c(1,1), 
lwd = c(1,1)), text = list(c("marginal mean profile",
"subject-specific profile")), columns = 2))

###Figure 6.6 
library(lattice)
trellis.device(color=F)
plot(model3.fit, resid(., type="p") ~ fitted(.) | factor(sicdegp), layout=c(3,1), aspect=2, abline=0)

plot(model3.fit, resid(., type="p") ~ fitted(.) | factor(sicdegp), id = 0.05, layout=c(3,1), aspect=2, abline=0)

###Figure 6.7 
plot(model3.fit, resid(.) ~ age.2, abline=0)

###Figure 6.8
qqnorm(model3.fit, ~resid(.) | factor(sicdegp), layout=c(3,1), aspect = 2, id = 0.05) 

detach(autism)

 