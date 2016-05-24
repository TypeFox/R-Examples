###########################################################
# Computation of Sensitivity Indexes (SI)
# for the so-called Ishigami function
# via Polynomial Chaos Expansion (PCE)
###########################################################


# Load of necessary functions
library("polychaosbasics")

nlhs<-200 # the number of rows of the complete dataset is 20000
degree<-6
#############################################
# Build Legendre polynomial
#############################################
set.seed(42)
pce <- analyticsPolyLeg(nlhs, degree, 'ishigami')
# -----------------------------------------
# Default display
# -----------------------------------------
print(pce)

# -----------------------------------------
# More display
# -----------------------------------------
print(pce, all=TRUE)

# -----------------------------------------
# All the components in the returned object
# -----------------------------------------
getNames(pce)

#############################################
# PCESI calculation
#############################################
retour<- PCESI(pce)
# -----------------------------------------
# Default display
# -----------------------------------------
print(retour)

# -----------------------------------------
# More display
# -----------------------------------------
print(retour, all=TRUE)

# -----------------------------------------
# All the components in the returned object
# -----------------------------------------
getNames(retour)

# -----------------------------------------
# Individual Monomial SI
# -----------------------------------------
print(retour@IMSI)

# -----------------------------------------
# Regression coefficients
# -----------------------------------------
print(retour@coef)


#############################################
# Plot  Y.hat (x-axis) against Y.hat (x-axis)
# and regression line
#############################################
y.hat <-retour@y.hat
y.obs <- pce[, "Y"]
reg <- lm(y.hat ~ y.obs)
jpeg(file="Ishigami.jpg")
 plot(y.hat, y.obs,
      xlab="metamodel output", ylab="computer model output",
      main="Ishigami test", sub="Scatter plot and regression line")
lines(reg$fitted.values, y.obs)
 dev.off()

